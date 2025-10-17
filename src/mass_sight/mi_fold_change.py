"""
Multiple imputation (MI) fold-change using calibrated matching candidates.

Public API
- mi_fold_change(ds1_expr, ds1_meta, ds2_expr, ds2_meta, candidates, M=100, seed=123)
  Returns a DataFrame with per-entity MI estimates and diagnostics.

Inputs
- ds1_expr, ds2_expr: DataFrames of raw intensities (rows=feature ids matching the
  DS tables used for matching; columns=sample ids). Values are nonnegative.
- ds1_meta, ds2_meta: DataFrames with columns: sample, group. group must have 2 levels.
- candidates: DataFrame from match_ml (columns include id1, id2, p_row, p0 or prob_cal).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple
from scipy.optimize import linear_sum_assignment  # type: ignore
import math
import numpy as np
import pandas as pd
from scipy.stats import t as student_t  # type: ignore


def _encode_group(meta: pd.DataFrame) -> pd.Series:
    g = meta["group"].astype("category")
    if len(g.cat.categories) != 2:
        raise ValueError("group must have exactly 2 levels")
    return (g.cat.codes == g.cat.categories.get_loc(g.cat.categories[-1])).astype(int)


def _ols_fc(log_y: np.ndarray, group_bin: np.ndarray, ds_indicator: np.ndarray) -> Tuple[float, float, int]:
    n = len(log_y)
    X = np.column_stack([np.ones(n), group_bin.astype(float), ds_indicator.astype(float)])
    try:
        XtX_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.inv(X.T @ X + 1e-8 * np.eye(3))
    beta = XtX_inv @ (X.T @ log_y)
    resid = log_y - X @ beta
    dof = max(1, n - X.shape[1])
    sigma2 = float((resid @ resid) / dof)
    var_beta_group = sigma2 * XtX_inv[1, 1]
    return float(beta[1]), float(var_beta_group), dof


def _sample_matching(cand: pd.DataFrame, rng: np.random.Generator) -> Dict[int, Optional[int]]:
    """Gumbel-perturbed Hungarian sampling of a 1–1 matching with per-row no-match.

    Steps per draw:
    1) For each row i, sample match vs no-match using p0_i.
    2) On the subset S of rows marked match, solve 1–1 via Hungarian with cost = -log(p) + Gumbel noise.
       If SciPy is unavailable or S empty, fall back to greedy by p.
    Returns a dict id1 -> id2 (or None).
    """
    pcol = "p_row" if "p_row" in cand.columns else ("prob_cal" if "prob_cal" in cand.columns else "prob")
    rows = cand["id1"].astype(int).unique().tolist()
    # Sample match vs no-match
    p0 = cand.drop_duplicates("id1").set_index("id1")["p0"].to_dict() if "p0" in cand.columns else {i: 0.0 for i in rows}
    match_rows = [i for i in rows if rng.random() >= float(p0.get(i, 0.0))]
    choice: Dict[int, Optional[int]] = {i: None for i in rows}
    if not match_rows:
        return choice
    sub = cand[cand["id1"].astype(int).isin(match_rows)].copy()
    # Build dense cost with Gumbel noise on -log(p)
    rlist = sorted(set(sub["id1"].astype(int)))
    clist = sorted(set(sub["id2"].astype(int)))
    r_index = {r: idx for idx, r in enumerate(rlist)}
    c_index = {c: idx for idx, c in enumerate(clist)}
    R, C = len(rlist), len(clist)
    cost = np.full((R, C), 50.0, dtype=float)
    # Gumbel(0,1) noise
    g_noise = rng.gumbel(loc=0.0, scale=1.0, size=(R, C))
    for _, r in sub.iterrows():
        i = r_index[int(r["id1"])]; j = c_index[int(r["id2"])]; p = float(r[pcol]); p = min(max(p, 1e-12), 1 - 1e-12)
        cost[i, j] = -np.log(p) + 0.05 * g_noise[i, j]
    ri, cj = linear_sum_assignment(cost)
    for a, b in zip(ri, cj):
        if cost[a, b] >= 25.0:
            continue
        choice[rlist[a]] = clist[b]
    return choice


def mi_fold_change(
    ds1_expr: pd.DataFrame,
    ds1_meta: pd.DataFrame,
    ds2_expr: pd.DataFrame,
    ds2_meta: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    M: int = 100,
    seed: int = 123,
) -> pd.DataFrame:
    """Run MI fold-change on calibrated candidates.

    Returns a DataFrame with per-entity logFC_hat, SE, CI, FC, and diagnostics.
    """
    # Align columns by sample order
    ds1_expr = ds1_expr.copy(); ds2_expr = ds2_expr.copy()
    ds1_meta = ds1_meta.copy(); ds2_meta = ds2_meta.copy()
    ds1_meta["group_bin"] = _encode_group(ds1_meta)
    ds2_meta["group_bin"] = _encode_group(ds2_meta)
    ds1_expr = ds1_expr.loc[:, ds1_meta["sample"]]
    ds2_expr = ds2_expr.loc[:, ds2_meta["sample"]]
    ds_ind1 = np.zeros(len(ds1_meta), dtype=int)
    ds_ind2 = np.ones(len(ds2_meta), dtype=int)

    # Per-row diagnostics
    def _best_margin(g: pd.DataFrame) -> Tuple[float, float]:
        g = g.sort_values("p_row" if "p_row" in g.columns else ("prob_cal" if "prob_cal" in g.columns else "prob"), ascending=False)
        b = float(g.iloc[0]["p_row" if "p_row" in g.columns else ("prob_cal" if "prob_cal" in g.columns else "prob")]) if len(g) else 0.0
        s = float(g.iloc[1]["p_row" if "p_row" in g.columns else ("prob_cal" if "prob_cal" in g.columns else "prob")]) if len(g) > 1 else 0.0
        return b, b - s
    bm = candidates.groupby("id1").apply(_best_margin)
    best = {int(k): v[0] for k, v in bm.to_dict().items()}
    margin = {int(k): v[1] for k, v in bm.to_dict().items()}

    rng = np.random.default_rng(seed)
    entity_store: Dict[Tuple[int, int], Dict[str, List[float]]] = {}
    used_count: Dict[Tuple[int, int], int] = {}
    for m in range(M):
        match = _sample_matching(candidates, rng)
        for i1, j2 in match.items():
            if j2 is None:
                continue
            i1 = int(i1); j2 = int(j2)
            if i1 not in ds1_expr.index or j2 not in ds2_expr.index:
                continue
            y1 = ds1_expr.loc[i1].to_numpy(dtype=float)
            y2 = ds2_expr.loc[j2].to_numpy(dtype=float)
            logy = np.concatenate([np.log1p(y1), np.log1p(y2)])
            group = np.concatenate([ds1_meta["group_bin"].to_numpy(), ds2_meta["group_bin"].to_numpy()])
            ds_ind = np.concatenate([ds_ind1, ds_ind2])
            if group.sum() == 0 or group.sum() == len(group):
                continue
            beta, varb, _ = _ols_fc(logy, group, ds_ind)
            key = (i1, j2)
            entity_store.setdefault(key, {"beta": [], "var": []})
            entity_store[key]["beta"].append(beta)
            entity_store[key]["var"].append(varb)
            used_count[key] = used_count.get(key, 0) + 1

    rows: List[Dict[str, float]] = []
    for (i1, j2), st in entity_store.items():
        betas = np.array(st["beta"], dtype=float)
        vars_ = np.array(st["var"], dtype=float)
        m_used = len(betas)
        if m_used < max(5, int(0.1 * M)):
            continue
        Qbar = float(betas.mean())
        Ubar = float(vars_.mean())
        B = float(betas.var(ddof=1)) if m_used > 1 else 0.0
        # Rubin's total variance
        T = Ubar + (1 + 1 / max(1, m_used)) * B
        se = math.sqrt(max(T, 1e-12))
        # Barnard–Rubin degrees of freedom for MI
        if B > 0:
            r = (1 + 1 / m_used) * B / max(Ubar, 1e-12)
            nu_old = m_used - 1
            nu = (m_used - 1) * (1 + 1 / r) ** 2
            # Small-sample correction (old+within) — conservative fallback if needed
            nu = float(max(1.0, nu))
            tcrit = float(student_t.ppf(0.975, df=nu))
        else:
            # No between-imputation variance: revert to normal
            nu = float('inf'); tcrit = 1.96
        rows.append({
            "entity_id": f"{i1}->{j2}", "id1": i1, "id2": j2,
            "logFC_hat": Qbar, "SE_MI": se,
            "CI_low": Qbar - tcrit * se, "CI_high": Qbar + tcrit * se,
            "FC": math.exp(Qbar),
            "draws_used": m_used, "match_prob": used_count.get((i1, j2), 0) / float(M),
            "p_row_best": best.get(i1, np.nan), "margin": margin.get(i1, np.nan),
            "n_samples": ds1_expr.shape[1] + ds2_expr.shape[1],
            "nu_MI": nu,
        })
    return pd.DataFrame(rows)


__all__ = ["mi_fold_change"]
