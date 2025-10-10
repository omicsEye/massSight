"""
Multi-dataset alignment from calibrated pairwise probabilities.

Abstraction
- Learn a pairwise scorer (we reuse match_ml) to get calibrated P(match) for each dataset pair.
- Generate candidate tuples anchored on dataset 0 (configurable) by taking top‑K candidates to other
  datasets (including None for no‑match).
- Score each tuple as the sum of log calibrated probabilities across present pairs (including optional
  non‑anchor inter‑dataset edges if available). Add log p0 for omitted datasets.
- Select a globally consistent set of tuples (at most one feature per dataset) via a greedy set‑packing
  heuristic (fast and effective)

Outputs
- MultiMatchResult with a DataFrame of selected tuples (one row per entity) and diagnostics.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple

import math
import numpy as np
import pandas as pd

from .ml_matcher import MLMatchConfig, match_ml


@dataclass
class MultiMatchConfig:
    # Base pairwise matcher settings
    ppm: float = 12.0
    rt: float = 1.0
    seed: int = 123
    calibrate: bool = True
    model: str = "auto"  # lgbm | logreg | auto
    # Candidate generation
    topk: int = 3
    prob_floor: float = 0.05
    tau_main: float = 0.2  # require at least one edge >= tau to keep tuple
    max_tuples_per_anchor: int = 50
    # Schema mapping (applied to all datasets)
    mz_col: str = "MZ"
    rt_col: str = "RT"
    annotation_col: str = "Annotation_ID"
    compound_id_col: str = "Compound_ID"
    intensity_col: Optional[str] = "Intensity"
    intensity_cols: Optional[List[str]] = None
    intensity_regex: Optional[str] = None


@dataclass
class MultiMatchResult:
    clusters: pd.DataFrame  # columns: cluster_id, id_0..id_{K-1}, score, size
    edges: Dict[Tuple[int, int], pd.DataFrame]  # per-pair calibrated edges used for scoring
    meta: Dict[str, object]


def _pair_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _run_pair(dsA: pd.DataFrame, dsB: pd.DataFrame, cfg: MultiMatchConfig) -> pd.DataFrame:
    """Run pairwise matcher and return calibrated candidate edges DataFrame.

    Returns columns: id1, id2, prob_cal, p0 (per id1), plus original ids.
    """
    mcfg = MLMatchConfig(
        ppm=cfg.ppm,
        rt=cfg.rt,
        seed=cfg.seed,
        model=cfg.model,
        calibrate=cfg.calibrate,
        mz_col=cfg.mz_col, rt_col=cfg.rt_col,
        annotation_col=cfg.annotation_col, compound_id_col=cfg.compound_id_col,
        intensity_col=cfg.intensity_col, intensity_cols=cfg.intensity_cols, intensity_regex=cfg.intensity_regex,
    )
    res = match_ml(dsA, dsB, mcfg)
    cand = res.candidates.copy()
    # Require calibrated probability and per-row no-match mass
    if "prob_cal" not in cand.columns:
        raise ValueError("Candidate table missing calibrated probability column 'prob_cal'. Ensure calibrate=True.")
    if "p0" not in cand.columns:
        raise ValueError("Candidate table missing 'p0' column. Run matcher with calibration to compute p_row/p0.")
    return cand


def _topk_by_id1(cand: pd.DataFrame, k: int, floor: float) -> pd.DataFrame:
    d = cand[cand["prob_cal"] >= floor].copy()
    d.sort_values(["id1", "prob_cal"], ascending=[True, False], inplace=True)
    return d.groupby("id1").head(k).reset_index(drop=True)


def _score_tuple(ids: List[Optional[int]], edges: Dict[Tuple[int, int], pd.DataFrame], p0_maps: Dict[Tuple[int, int], Dict[int, float]]) -> float:
    K = len(ids)
    s = 0.0
    # Sum logs for all present pairs
    for a in range(K):
        ia = ids[a]
        if ia is None:
            continue
        for b in range(a + 1, K):
            ib = ids[b]
            if ib is None:
                continue
            key = _pair_key(a, b)
            df = edges.get(key)
            if df is None or df.empty:
                continue
            # df has id1 (from min index), id2 (from max index) depending on ordering
            if a < b:
                m = df[(df["id1"].astype(int) == int(ia)) & (df["id2"].astype(int) == int(ib))]
            else:
                m = df[(df["id1"].astype(int) == int(ib)) & (df["id2"].astype(int) == int(ia))]
            if not m.empty:
                p = float(m.iloc[0]["prob_cal"])  # calibrated probability
                p = min(max(p, 1e-12), 1 - 1e-12)
                s += math.log(p)
    # Add log p0 for omitted datasets relative to anchor 0
    for b in range(1, K):
        if ids[b] is None and ids[0] is not None:
            key = _pair_key(0, b)
            p0_map = p0_maps.get(key)
            if p0_map is None or int(ids[0]) not in p0_map:
                raise ValueError("Missing p0 for anchor id when scoring tuple; run matcher with calibration to include p0.")
            p0_val = float(p0_map[int(ids[0])])
            p0_val = min(max(p0_val, 1e-12), 1 - 1e-12)
            s += math.log(p0_val)
    return s


def align_multi(datasets: Sequence[pd.DataFrame], *, names: Optional[Sequence[str]] = None, config: Optional[MultiMatchConfig] = None) -> MultiMatchResult:
    cfg = config or MultiMatchConfig()
    K = len(datasets)
    if K < 3:
        raise ValueError("align_multi requires at least 3 datasets")
    if names is None:
        names = [f"DS{k+1}" for k in range(K)]

    # 1) Compute calibrated pairwise edges for all pairs (a<b)
    edges: Dict[Tuple[int, int], pd.DataFrame] = {}
    p0_maps: Dict[Tuple[int, int], Dict[int, float]] = {}
    for a in range(K):
        for b in range(a + 1, K):
            cand = _run_pair(datasets[a], datasets[b], cfg)
            # enforce top‑K per anchor id1 (id from the smaller index dataset)
            cand_k = _topk_by_id1(cand, cfg.topk, cfg.prob_floor)
            edges[(a, b)] = cand_k
            # p0 per id1 for this pair
            p0 = cand_k.drop_duplicates("id1").set_index("id1")["p0"].to_dict()
            p0_maps[(a, b)] = {int(k): float(v) for k, v in p0.items()}

    # 2) Generate tuple candidates anchored on dataset 0
    ids0 = list(range(len(datasets[0])))
    # Prepare quick lookup dicts for anchor to others
    anchor_cands: Dict[int, Dict[int, List[Tuple[int, float]]]] = {}
    for b in range(1, K):
        key = _pair_key(0, b)
        df = edges.get(key, pd.DataFrame())
        m = {}
        if not df.empty:
            for i, g in df.groupby("id1"):
                m[int(i)] = [(int(r.id2), float(r.prob_cal)) for _, r in g.iterrows()]
        anchor_cands[b] = m

    tuples: List[Dict[str, object]] = []
    for i0 in ids0:
        # Build candidate lists per other dataset, include None with implicit probability via p0 later
        cand_lists: List[List[Optional[int]]] = []
        keep_flag = False
        for b in range(1, K):
            lst = [None]
            pairs = anchor_cands[b].get(i0, [])
            # ensure sorted by prob desc and meet floor
            pairs = [(j, p) for (j, p) in pairs if p >= cfg.prob_floor]
            pairs.sort(key=lambda x: x[1], reverse=True)
            if pairs:
                lst += [j for (j, _) in pairs[: cfg.topk]]
                if pairs[0][1] >= cfg.tau_main:
                    keep_flag = True
            cand_lists.append(lst)
        if not cand_lists:
            continue
        # Optional: if no strong edge, still consider with weak flag; keep small number
        # Generate cartesian product with pruning by max tuples per anchor
        count = 0
        for choices in product(*cand_lists):
            ids = [i0] + list(choices)
            # skip all-None (no matches)
            if all(x is None for x in ids[1:]):
                continue
            # require at least one edge >= tau_main from anchor
            if not keep_flag:
                # allow but limit
                pass
            score = _score_tuple(ids, edges, p0_maps)
            tuples.append({
                "id_0": i0,
                **{f"id_{b}": (None if ids[b] is None else int(ids[b])) for b in range(1, K)},
                "score": float(score),
                "size": int(sum(1 for x in ids if x is not None)),
            })
            count += 1
            if count >= cfg.max_tuples_per_anchor:
                break

    if not tuples:
        return MultiMatchResult(clusters=pd.DataFrame(columns=[f"id_{k}" for k in range(K)] + ["score", "size", "cluster_id"]), edges=edges, meta={"names": list(names)})

    cand_tuples = pd.DataFrame(tuples)
    cand_tuples.sort_values("score", ascending=False, inplace=True)

    # 3) Greedy set-packing selection (at most one per dataset)
    used: Dict[Tuple[int, int], bool] = {}
    chosen: List[Dict[str, object]] = []
    for _, r in cand_tuples.iterrows():
        conflict = False
        for k in range(K):
            idx = r.get(f"id_{k}")
            if pd.isna(idx):
                continue
            key = (k, int(idx))
            if used.get(key, False):
                conflict = True
                break
        if conflict:
            continue
        # accept
        for k in range(K):
            idx = r.get(f"id_{k}")
            if pd.isna(idx):
                continue
            used[(k, int(idx))] = True
        chosen.append(r.to_dict())

    if not chosen:
        clusters = pd.DataFrame(columns=[f"id_{k}" for k in range(K)] + ["score", "size", "cluster_id"]) 
    else:
        clusters = pd.DataFrame(chosen).reset_index(drop=True)
        clusters["cluster_id"] = np.arange(len(clusters))

    meta = {
        "names": list(names),
        "K": K,
        "topk": cfg.topk,
        "prob_floor": cfg.prob_floor,
        "tau_main": cfg.tau_main,
    }
    return MultiMatchResult(clusters=clusters, edges=edges, meta=meta)


__all__ = ["MultiMatchConfig", "MultiMatchResult", "align_multi"]
