from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .lcms_utils import ADDUCT_MASS_NEG, ADDUCT_MASS_POS, PROTON_MASS


@dataclass(frozen=True)
class AdductGroupingConfig:
    polarity: str  # "positive" | "negative"
    adducts: str = "expanded"  # "conservative" | "expanded"
    delta_strategy: str = "base_to_others"  # "base_to_others" | "all_pairs"
    include_neutral_mass: bool = False
    ppm: float = 5.0
    rt_window_min: float = 0.1
    include_isotope: bool = False
    isotope_shift_da: float = 1.0033548378
    max_edges_per_feature: int = 200


def _normalize_polarity(polarity: str) -> str:
    pol = str(polarity or "").strip().lower()
    if pol in {"pos", "positive"}:
        return "positive"
    if pol in {"neg", "negative"}:
        return "negative"
    raise ValueError(f"Unsupported polarity: {polarity!r}")


def _adduct_masses(cfg: AdductGroupingConfig) -> List[Tuple[str, float]]:
    pol = _normalize_polarity(cfg.polarity)
    mode = str(cfg.adducts or "expanded").strip().lower()
    if mode not in {"conservative", "expanded"}:
        raise ValueError(f"Unsupported adducts mode: {cfg.adducts!r}")

    out: List[Tuple[str, float]] = []
    if bool(getattr(cfg, "include_neutral_mass", False)):
        out.append(("NEUTRAL", 0.0))
    if pol == "positive":
        # Observed m/z ~= neutral + adduct_mass
        if mode == "conservative":
            out.append(("H", float(PROTON_MASS)))
        else:
            out.extend([(k, float(v)) for k, v in ADDUCT_MASS_POS.items()])
    else:
        # Represent [M-H]- as an adduct mass of -H (so neutral = mz - (-H) = mz + H).
        out.append(("M_H", -float(PROTON_MASS)))
        if mode == "expanded":
            out.extend([(k, float(v)) for k, v in ADDUCT_MASS_NEG.items()])
    return out


def _adduct_deltas(cfg: AdductGroupingConfig) -> List[Tuple[float, Tuple[str, str]]]:
    masses = _adduct_masses(cfg)
    if not masses:
        raise ValueError("No adduct masses configured for within-study grouping.")
    strat = str(getattr(cfg, "delta_strategy", "base_to_others") or "base_to_others").strip().lower()
    if strat not in {"base_to_others", "all_pairs"}:
        raise ValueError(f"Unsupported delta_strategy: {cfg.delta_strategy!r}")
    pol = _normalize_polarity(cfg.polarity)
    base_name = "H" if pol == "positive" else "M_H"
    base = [x for x in masses if x[0] == base_name]
    if not base:
        base = [masses[0]]
    base_mass = float(base[0][1])

    deltas: List[Tuple[float, Tuple[str, str]]] = []
    if strat == "all_pairs":
        for name1, m1 in masses:
            for name2, m2 in masses:
                d = float(m2 - m1)
                if not np.isfinite(d) or abs(d) <= 1e-12:
                    continue
                deltas.append((d, (name1, name2)))
    else:
        for name2, m2 in masses:
            if name2 == base_name:
                continue
            d = float(m2 - base_mass)
            if np.isfinite(d) and abs(d) > 1e-12:
                deltas.append((d, (base_name, name2)))
                deltas.append((-d, (name2, base_name)))
    if bool(cfg.include_isotope):
        s = float(cfg.isotope_shift_da)
        if np.isfinite(s) and s > 0:
            deltas.append((s, ("M", "M_C13")))
            deltas.append((-s, ("M_C13", "M")))
    # Dedupe within numerical tolerance.
    deltas.sort(key=lambda x: x[0])
    out: List[Tuple[float, Tuple[str, str]]] = []
    for d, lab in deltas:
        if out and abs(d - out[-1][0]) <= 1e-9:
            continue
        out.append((d, lab))
    return out


def build_adduct_edges(
    mz: np.ndarray,
    rt: np.ndarray,
    *,
    cfg: AdductGroupingConfig,
) -> pd.DataFrame:
    """
    Build within-study candidate adduct/isotope edges using:
      - co-elution (|ΔRT| <= rt_window_min)
      - m/z differences matching an adduct mass delta (within ppm)

    Returns an edge table with columns:
      id1, id2, delta_da, ppm_error, rt_diff_min, adduct_from, adduct_to
    """
    mz = np.asarray(mz, dtype=float)
    rt = np.asarray(rt, dtype=float)
    if mz.ndim != 1 or rt.ndim != 1 or mz.size != rt.size:
        raise ValueError("mz and rt must be 1D aligned arrays.")

    ppm = float(cfg.ppm)
    if not np.isfinite(ppm) or ppm <= 0:
        raise ValueError("cfg.ppm must be finite and > 0.")
    rt_w = float(cfg.rt_window_min)
    if not np.isfinite(rt_w) or rt_w <= 0:
        raise ValueError("cfg.rt_window_min must be finite and > 0.")

    n = int(mz.size)
    finite = np.isfinite(mz) & (mz > 0) & np.isfinite(rt)
    order = np.argsort(np.where(finite, mz, np.inf), kind="mergesort")
    mz_s = mz[order]

    deltas = _adduct_deltas(cfg)
    max_edges = int(cfg.max_edges_per_feature)
    if max_edges <= 0:
        max_edges = 1

    # Keep the best (lowest ppm error) match per unordered pair.
    best: Dict[Tuple[int, int], Tuple[float, int, int, float, float, str, str]] = {}
    # key -> (ppm_error, id1, id2, delta_da, rt_diff_min, adduct_from, adduct_to)

    for i in range(n):
        if not finite[i]:
            continue
        mi = float(mz[i])
        rti = float(rt[i])
        edge_count = 0
        for d, (a1, a2) in deltas:
            target = mi + float(d)
            if not np.isfinite(target) or target <= 0:
                continue
            tol = (ppm * 1e-6) * target
            lo = target - tol
            hi = target + tol
            l = int(np.searchsorted(mz_s, lo, side="left"))
            r = int(np.searchsorted(mz_s, hi, side="right"))
            if r <= l:
                continue
            cand = order[l:r]
            for j in cand.tolist():
                if j == i or not finite[j]:
                    continue
                rtj = float(rt[j])
                rt_diff = abs(rtj - rti)
                if rt_diff > rt_w:
                    continue
                mj = float(mz[j])
                ppm_err = abs(mj - target) / max(target, 1.0) * 1e6
                if ppm_err > ppm:
                    continue
                # Unordered key for dedupe.
                a, b = (i, j) if i < j else (j, i)
                key = (a, b)
                prev = best.get(key)
                if prev is None or float(ppm_err) < float(prev[0]):
                    best[key] = (float(ppm_err), int(i), int(j), float(d), float(rt_diff), str(a1), str(a2))
                edge_count += 1
                if edge_count >= max_edges:
                    break
            if edge_count >= max_edges:
                break

    rows = list(best.values())
    if not rows:
        return pd.DataFrame(
            columns=[
                "id1",
                "id2",
                "delta_da",
                "ppm_error",
                "rt_diff_min",
                "adduct_from",
                "adduct_to",
            ]
        )
    df = pd.DataFrame(
        rows,
        columns=[
            "ppm_error",
            "id1",
            "id2",
            "delta_da",
            "rt_diff_min",
            "adduct_from",
            "adduct_to",
        ],
    )
    df = df.loc[:, ["id1", "id2", "delta_da", "ppm_error", "rt_diff_min", "adduct_from", "adduct_to"]]
    df = df.sort_values(["id1", "id2"], kind="mergesort").reset_index(drop=True)
    return df


def _union_find(n: int, edges: np.ndarray) -> np.ndarray:
    parent = np.arange(n, dtype=int)
    rank = np.zeros(n, dtype=int)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for a, b in edges.tolist():
        union(int(a), int(b))

    roots = np.array([find(i) for i in range(n)], dtype=int)
    # Compress to 0..G-1.
    _, inv = np.unique(roots, return_inverse=True)
    return inv.astype(int)


def assign_adduct_groups(
    mz: np.ndarray,
    rt: np.ndarray,
    *,
    cfg: AdductGroupingConfig,
) -> Tuple[np.ndarray, pd.DataFrame, dict]:
    edges = build_adduct_edges(mz, rt, cfg=cfg)
    n = int(np.asarray(mz).size)
    if edges.empty:
        return np.arange(n, dtype=int), edges, {"n_edges": 0, "n_groups": n}

    pairs = edges.loc[:, ["id1", "id2"]].to_numpy(dtype=int, copy=False)
    group_id = _union_find(n, pairs)

    diag = {
        "n_features": int(n),
        "n_edges": int(edges.shape[0]),
        "n_groups": int(np.unique(group_id).size),
        "cfg": {
            "polarity": _normalize_polarity(cfg.polarity),
            "adducts": str(cfg.adducts),
            "delta_strategy": str(getattr(cfg, "delta_strategy", "base_to_others")),
            "include_neutral_mass": bool(getattr(cfg, "include_neutral_mass", False)),
            "ppm": float(cfg.ppm),
            "rt_window_min": float(cfg.rt_window_min),
            "include_isotope": bool(cfg.include_isotope),
            "max_edges_per_feature": int(cfg.max_edges_per_feature),
        },
    }
    return group_id, edges, diag


def infer_group_neutral_mass_and_adducts(
    mz: np.ndarray,
    group_id: np.ndarray,
    *,
    cfg: AdductGroupingConfig,
    ppm_tol: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Infer a neutral mass and adduct label per feature within each adduct group.

    Returns:
      neutral_mass_est (float, NaN if not assigned),
      adduct_label (object),
      neutral_ppm_error (float, NaN if not assigned)
    """
    mz = np.asarray(mz, dtype=float)
    gid = np.asarray(group_id, dtype=int)
    if mz.ndim != 1 or gid.ndim != 1 or mz.size != gid.size:
        raise ValueError("mz and group_id must be aligned 1D arrays.")

    ppm_tol = float(ppm_tol)
    if not np.isfinite(ppm_tol) or ppm_tol <= 0:
        raise ValueError("ppm_tol must be finite and > 0.")

    adducts = _adduct_masses(cfg)
    names = [a[0] for a in adducts]
    masses = np.array([a[1] for a in adducts], dtype=float)
    n = int(mz.size)

    neutral_est = np.full(n, np.nan, dtype=float)
    neutral_err_ppm = np.full(n, np.nan, dtype=float)
    labels = np.full(n, "", dtype=object)

    for g in np.unique(gid).tolist():
        idx = np.flatnonzero(gid == int(g))
        if idx.size <= 1:
            # Singleton: no labeling attempt.
            continue
        mz_g = mz[idx]
        if not np.all(np.isfinite(mz_g)) or np.any(mz_g <= 0):
            continue

        # Candidate neutrals from all feature×adduct combinations.
        neutral_cand = (mz_g[:, None] - masses[None, :]).reshape(-1)
        # Score each candidate by how many distinct features can be assigned.
        best_score = -1
        best_neutral = float("nan")
        best_med_ppm = float("inf")

        for nm in neutral_cand.tolist():
            nm = float(nm)
            if (not np.isfinite(nm)) or nm <= 0:
                continue
            # For each feature, see if any adduct matches within ppm.
            ppm_err = np.abs((mz_g[:, None] - masses[None, :]) - nm) / nm * 1e6
            min_ppm = np.min(ppm_err, axis=1)
            ok = min_ppm <= ppm_tol
            score = int(np.sum(ok))
            if score < best_score:
                continue
            med_ppm = float(np.median(min_ppm[ok])) if score > 0 else float("inf")
            if score > best_score or (score == best_score and med_ppm < best_med_ppm):
                best_score = score
                best_neutral = nm
                best_med_ppm = med_ppm

        if best_score <= 0 or (not np.isfinite(best_neutral)) or best_neutral <= 0:
            continue

        # Assign each feature an adduct label based on best_neutral.
        ppm_err = np.abs((mz_g[:, None] - masses[None, :]) - best_neutral) / best_neutral * 1e6
        best_adduct = np.argmin(ppm_err, axis=1)
        best_ppm = ppm_err[np.arange(idx.size), best_adduct]
        ok = best_ppm <= ppm_tol
        if not np.any(ok):
            continue

        # Re-estimate neutral mass from assigned members for robustness.
        neutral_vals = mz_g[ok] - masses[best_adduct[ok]]
        neutral_hat = float(np.median(neutral_vals))
        if not np.isfinite(neutral_hat) or neutral_hat <= 0:
            continue

        # Store per-feature assignments.
        for pos, fi in enumerate(idx.tolist()):
            if not ok[pos]:
                continue
            ai = int(best_adduct[pos])
            labels[fi] = str(names[ai])
            nm_i = float(mz[fi] - float(masses[ai]))
            neutral_est[fi] = float(neutral_hat)
            neutral_err_ppm[fi] = float(abs(nm_i - neutral_hat) / neutral_hat * 1e6)

    return neutral_est, labels, neutral_err_ppm
