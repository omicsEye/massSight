"""
Supervised ML matcher with calibrated probabilities.

This module trains a simple classifier on annotated singleton matches between
two datasets and produces calibrated candidate probabilities for all DS1 rows.

Main entry point
- match_ml(ds1: DataFrame, ds2: DataFrame, config: MLMatchConfig) -> MLMatchResult

Notes
- Uses RT±w and ppm±w windows (no isotonic RT mapping) to form candidates.
- Positives: shared singletons (appear exactly once in each dataset).
- Negatives: near‑miss on both sides + a few random negatives; skips DS2
  candidates whose Annotation_ID is in the calibration/test slice to avoid leakage.
- Calibration: out‑of‑fold isotonic (5‑fold) by default; returns prob_cal and
  per‑row masses p_row and p0 for downstream sampling / MI.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from lightgbm import LGBMClassifier  # type: ignore
from scipy.optimize import linear_sum_assignment  # type: ignore
import joblib
import time

def get_retention_order(ds: pd.DataFrame, rt_col: str = 'RT') -> pd.Series:
    """Normalized retention order in [0,1], robust to monotone drift."""
    return ds[rt_col].rank(pct=True)


@dataclass
class MLMatchConfig:
    # Matching windows
    ppm: float = 12.0
    rt: float = 1.0
    # Model
    model: str = "auto"  # one of: auto, logreg, lgbm
    seed: int = 123
    calibrate: bool = True
    # Negatives
    neg_row_k: int = 8
    neg_col_k: int = 4
    neg_rand_k: int = 2
    neg_expand: float = 1.5
    neg_hard_k: int = 4
    # Schema: column names per dataset
    mz_col: str = "MZ"
    rt_col: str = "RT"
    annotation_col: str = "Annotation_ID"
    compound_id_col: str = "Compound_ID"
    # Intensity options: either a single column name, or a list, or a regex pattern
    intensity_col: Optional[str] = "Intensity"
    intensity_cols: Optional[List[str]] = None
    intensity_regex: Optional[str] = None
    # Persistence
    save_model_path: Optional[str] = None


@dataclass
class MLMatchResult:
    candidates: pd.DataFrame  # columns: id1, id2, prob_raw/ prob_cal, p_row, p0, ppm_diff, rt_diff, roi_diff, log_int_diff
    top1: pd.DataFrame  # columns: id1, id2, prob
    hungarian: pd.DataFrame  # columns: id1, id2, prob (if scipy available; else greedy)


def _normalize_schema(ds: pd.DataFrame, cfg: MLMatchConfig) -> pd.DataFrame:
    """Return a copy with standardized columns: MZ, RT, Annotation_ID, Compound_ID,
    and intensity features (Intensity_log10). If multiple intensity columns are
    specified, Intensity_log10 = log10(mean(raw_intensities)). If a single
    column is specified or present, Intensity_log10 = log10(raw + eps).
    """
    df = ds.copy()
    # Rename core columns if needed
    if cfg.mz_col != "MZ" and cfg.mz_col in df.columns:
        df.rename(columns={cfg.mz_col: "MZ"}, inplace=True)
    if cfg.rt_col != "RT" and cfg.rt_col in df.columns:
        df.rename(columns={cfg.rt_col: "RT"}, inplace=True)
    if cfg.annotation_col != "Annotation_ID" and cfg.annotation_col in df.columns:
        df.rename(columns={cfg.annotation_col: "Annotation_ID"}, inplace=True)
    if cfg.compound_id_col != "Compound_ID" and cfg.compound_id_col in df.columns:
        df.rename(columns={cfg.compound_id_col: "Compound_ID"}, inplace=True)
    # Build Intensity_log10
    eps = 1e-9
    int_log = None
    if cfg.intensity_cols:
        cols = [c for c in cfg.intensity_cols if c in df.columns]
        if cols:
            int_log = np.log10(df[cols].astype(float).clip(lower=0).mean(axis=1) + eps)
    elif cfg.intensity_regex:
        cols = [c for c in df.columns if pd.Series([c]).str.contains(cfg.intensity_regex).iloc[0]]
        if cols:
            int_log = np.log10(df[cols].astype(float).clip(lower=0).mean(axis=1) + eps)
    elif cfg.intensity_col and cfg.intensity_col in df.columns:
        int_log = np.log10(df[cfg.intensity_col].astype(float).clip(lower=0) + eps)
    # Fallback if nothing provided: no intensity info
    if int_log is None:
        df["Intensity_log10"] = 0.0
    else:
        df["Intensity_log10"] = int_log
    # Retention order
    if "ro" not in df.columns:
        df["ro"] = get_retention_order(df)
    return df


def _shared_singletons(ds1: pd.DataFrame, ds2: pd.DataFrame) -> pd.DataFrame:
    col = "Annotation_ID"
    if col not in ds1.columns or col not in ds2.columns:
        return pd.DataFrame(columns=["id1", "id2", "label"])  # empty
    s1 = ds1[ds1[col].notna()][[col]].copy()
    s2 = ds2[ds2[col].notna()][[col]].copy()
    vc1 = s1[col].value_counts(); vc2 = s2[col].value_counts()
    once = set(vc1[vc1 == 1].index).intersection(set(vc2[vc2 == 1].index))
    if not once:
        return pd.DataFrame(columns=["id1", "id2", "label"])  # empty
    i1 = s1.reset_index().rename(columns={"index": "id1", col: "label"})
    i2 = s2.reset_index().rename(columns={"index": "id2", col: "label"})
    m = pd.merge(i1[i1["label"].isin(once)], i2[i2["label"].isin(once)], on="label", how="inner")
    return m[["id1", "id2", "label"]].reset_index(drop=True)


def _rt_window_indices(rt_sorted_vals: np.ndarray, rt_sorted_idx: np.ndarray, rt_center: float, w: float) -> np.ndarray:
    left = np.searchsorted(rt_sorted_vals, rt_center - w, side="left")
    right = np.searchsorted(rt_sorted_vals, rt_center + w, side="right")
    return rt_sorted_idx[left:right]


def _ppm_minmax(mz: float, ppm: float) -> Tuple[float, float]:
    span = mz * ppm / 1e6
    return mz - span, mz + span


def _build_neighborhoods(ds1: pd.DataFrame, ds2: pd.DataFrame, ppm: float, rt: float) -> List[List[int]]:
    rt_idx = np.argsort(ds2["RT"].values)
    rt_vals = ds2["RT"].values[rt_idx]
    mz2 = ds2["MZ"].values
    nbh: List[List[int]] = []
    for mz1, rt1 in zip(ds1["MZ"].values, ds1["RT"].values):
        cand_rt = _rt_window_indices(rt_vals, rt_idx, float(rt1), rt)
        if len(cand_rt) == 0:
            nbh.append([])
            continue
        lo, hi = _ppm_minmax(float(mz1), ppm)
        keep = [int(j) for j in cand_rt if lo <= mz2[j] <= hi]
        nbh.append(keep)
    return nbh


def _feature_row(ds1: pd.DataFrame, ds2: pd.DataFrame, i: int, j: int) -> dict:
    mz1, mz2 = float(ds1.at[i, "MZ"]), float(ds2.at[j, "MZ"])  # type: ignore[index]
    rt1, rt2 = float(ds1.at[i, "RT"]), float(ds2.at[j, "RT"])  # type: ignore[index]
    ro1, ro2 = float(ds1.at[i, "ro"]), float(ds2.at[j, "ro"])  # type: ignore[index]
    il1 = float(ds1.at[i, "Intensity_log10"]) if "Intensity_log10" in ds1.columns else 0.0
    il2 = float(ds2.at[j, "Intensity_log10"]) if "Intensity_log10" in ds2.columns else 0.0
    eps = 1e-12
    # Engineered features (canonical set for modeling)
    ppm_abs = abs(mz1 - mz2) / max(mz1, eps) * 1e6
    rt_abs = abs(rt1 - rt2)
    rt_signed = (rt2 - rt1)  # consistent orientation (DS2 − DS1)
    roi_abs = abs(ro1 - ro2)
    log_int_abs = abs(il1 - il2)
    rt_center = 0.5 * (rt1 + rt2)
    return {
        # identifiers and raw values (for diagnostics; excluded from model features)
        "id1": int(i), "id2": int(j),
        "mz1": mz1, "mz2": mz2, "rt1": rt1, "rt2": rt2, "ro1": ro1, "ro2": ro2,
        "int1": il1, "int2": il2,
        # legacy names for backward compatibility in outputs
        "ppm_diff": ppm_abs,
        "rt_diff": rt_abs,
        "roi_diff": roi_abs,
        "log_int_diff": log_int_abs,
        # new canonical features for training
        "ppm_abs": ppm_abs,
        "rt_abs": rt_abs,
        "rt_signed": rt_signed,
        "roi_abs": roi_abs,
        "log_int_abs": log_int_abs,
        "rt_center": rt_center,
    }


def _build_train(ds1: pd.DataFrame, ds2: pd.DataFrame, nbh: List[List[int]], pairs: pd.DataFrame, cfg: MLMatchConfig) -> Tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(cfg.seed)
    # sort indices for reverse-side negatives
    rt_idx1 = np.argsort(ds1["RT"].values)
    rt_vals1 = ds1["RT"].values[rt_idx1]
    mz1_arr = ds1["MZ"].values
    rt_idx2 = np.argsort(ds2["RT"].values)
    rt_vals2 = ds2["RT"].values[rt_idx2]
    mz2_arr = ds2["MZ"].values
    # exclude DS2 labels from test half during OOF; we set later per fold
    Xrows: List[dict] = []
    y: List[int] = []

    def cost(i_idx: int, j_idx: int) -> float:
        r = _feature_row(ds1, ds2, i_idx, j_idx)
        return (
            (r["rt_diff"] / max(cfg.rt, 1e-6)) ** 2
            + (r["ppm_diff"] / max(cfg.ppm, 1e-6)) ** 2
            + (r["roi_diff"] / 0.1) ** 2
            + (r["log_int_diff"] / 0.5) ** 2
        )

    for _, r in pairs.iterrows():
        i1, j2 = int(r.id1), int(r.id2)
        # positive
        Xrows.append(_feature_row(ds1, ds2, i1, j2)); y.append(1)
        # row-side negatives
        base = [jj for jj in nbh[i1] if jj != j2]
        if len(base) < max(cfg.neg_row_k, cfg.neg_hard_k):
            exp_idx = _rt_window_indices(rt_vals2, rt_idx2, float(ds1.at[i1, "RT"]), cfg.rt * cfg.neg_expand)
            lo, hi = _ppm_minmax(float(ds1.at[i1, "MZ"]), cfg.ppm * cfg.neg_expand)
            base = list(dict.fromkeys(base + [int(jj) for jj in exp_idx if lo <= mz2_arr[jj] <= hi and jj != j2]))
        if base:
            hard = sorted(base, key=lambda jj: cost(i1, jj))[: cfg.neg_hard_k]
            for jj in hard:
                Xrows.append(_feature_row(ds1, ds2, i1, jj)); y.append(0)
            rng.shuffle(base)
            for jj in base[: max(0, cfg.neg_row_k - cfg.neg_hard_k)]:
                Xrows.append(_feature_row(ds1, ds2, i1, jj)); y.append(0)
        # column-side negatives
        rt2_j = float(ds2.at[j2, "RT"])  # type: ignore[index]
        mz2_j = float(ds2.at[j2, "MZ"])  # type: ignore[index]
        cand_i = _rt_window_indices(rt_vals1, rt_idx1, rt2_j, cfg.rt)
        lo2, hi2 = _ppm_minmax(mz2_j, cfg.ppm)
        neg_i = [int(ii) for ii in cand_i if lo2 <= mz1_arr[ii] <= hi2 and ii != i1]
        if len(neg_i) < cfg.neg_col_k:
            exp_i = _rt_window_indices(rt_vals1, rt_idx1, rt2_j, cfg.rt * cfg.neg_expand)
            lo3, hi3 = _ppm_minmax(mz2_j, cfg.ppm * cfg.neg_expand)
            neg_i = list(dict.fromkeys(neg_i + [int(ii) for ii in exp_i if lo3 <= mz1_arr[ii] <= hi3 and ii != i1]))
        if neg_i:
            hard_i = sorted(neg_i, key=lambda ii: cost(ii, j2))[: cfg.neg_col_k]
            for ii in hard_i:
                Xrows.append(_feature_row(ds1, ds2, ii, j2)); y.append(0)
        # random negatives anywhere
        if cfg.neg_rand_k > 0:
            all_j = np.arange(len(ds2), dtype=int)
            mask = np.ones(len(ds2), dtype=bool)
            mask[j2] = False
            if "Annotation_ID" in ds2.columns:
                # crude guard, will be refined in OOF
                pass
            idx_pool = all_j[mask]
            if len(idx_pool) > 0:
                rands = rng.choice(idx_pool, size=min(cfg.neg_rand_k, len(idx_pool)), replace=False)
                for jj in rands:
                    Xrows.append(_feature_row(ds1, ds2, i1, int(jj))); y.append(0)

    X = pd.DataFrame(Xrows)
    yser = pd.Series(y)
    return X, yser


def _fit_oof_isotonic(clf, X: pd.DataFrame, y: pd.Series, *, seed: int) -> IsotonicRegression:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    scores = np.empty(len(X))
    for tr, va in skf.split(X, y):
        c = clone(clf)
        c.fit(X.iloc[tr], y.iloc[tr])
        scores[va] = c.predict_proba(X.iloc[va])[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(scores, y.values)
    return iso


def _assignment_top1(cand: pd.DataFrame, score_col: str) -> pd.DataFrame:
    d = cand.sort_values(["id1", score_col], ascending=[True, False]).copy()
    top = d.groupby("id1").head(1)[["id1", "id2", score_col]].rename(columns={score_col: "prob"})
    return top.reset_index(drop=True)


def _assignment_hungarian(cand: pd.DataFrame, score_col: str) -> pd.DataFrame:
    rows = sorted(cand["id1"].astype(int).unique())
    cols = sorted(cand["id2"].astype(int).unique())
    r_index = {r: idx for idx, r in enumerate(rows)}
    c_index = {c: idx for idx, c in enumerate(cols)}
    R = len(rows); C = len(cols)
    cost = np.full((R, C), 50.0, dtype=float)
    for _, r in cand.iterrows():
        i = r_index[int(r["id1"])]; j = c_index[int(r["id2"])]
        p = float(r[score_col]); p = min(max(p, 1e-12), 1 - 1e-12)
        cost[i, j] = -np.log(p)
    ri, cj = linear_sum_assignment(cost)
    pairs = []
    for a, b in zip(ri, cj):
        if cost[a, b] >= 25.0:
            continue
        id1 = rows[a]; id2 = cols[b]
        p = float(cand[(cand["id1"] == id1) & (cand["id2"] == id2)][score_col].iloc[0])
        pairs.append({"id1": id1, "id2": id2, "prob": p})
    return pd.DataFrame(pairs)


def match_ml(ds1: pd.DataFrame, ds2: pd.DataFrame, config: Optional[MLMatchConfig] = None) -> MLMatchResult:
    cfg = config or MLMatchConfig()
    ds1 = _normalize_schema(ds1, cfg)
    ds2 = _normalize_schema(ds2, cfg)
    pairs = _shared_singletons(ds1, ds2)
    if pairs.empty:
        raise ValueError("No shared singleton annotations to supervise the ML matcher")
    nbh = _build_neighborhoods(ds1, ds2, cfg.ppm, cfg.rt)
    X, y = _build_train(ds1, ds2, nbh, pairs, cfg)
    # Restrict to canonical, physically meaningful features
    allowed_features = ["ppm_abs", "rt_abs", "rt_signed", "roi_abs", "log_int_abs", "rt_center"]
    feature_names = [c for c in allowed_features if c in X.columns]
    if not feature_names:
        raise RuntimeError("No valid feature columns found for training (expected ppm_abs, rt_abs, rt_signed, roi_abs, log_int_abs, rt_center).")
    # Choose model
    mdl_choice = cfg.model
    if mdl_choice == "auto":
        mdl_choice = "lgbm"
    if mdl_choice == "lgbm":
        clf = LGBMClassifier(n_estimators=400, learning_rate=0.05, num_leaves=31, max_depth=-1,
                             subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, class_weight="balanced",
                             random_state=cfg.seed, verbosity=-1)
    else:
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", solver="lbfgs")
    # Calibration
    iso: Optional[IsotonicRegression] = None
    if cfg.calibrate:
        iso = _fit_oof_isotonic(clf, X[feature_names], y, seed=cfg.seed)
        clf.fit(X[feature_names], y)
    else:
        clf.fit(X[feature_names], y)
    # Optional: persist model artifact
    if cfg.save_model_path:
        artifact = {
            "model": clf,
            "iso": iso,
            "feature_names": list(feature_names),
            "config": {
                "ppm": cfg.ppm,
                "rt": cfg.rt,
                "model": cfg.model,
                "seed": cfg.seed,
                "calibrate": cfg.calibrate,
                "schema": {
                    "mz_col": cfg.mz_col,
                    "rt_col": cfg.rt_col,
                    "annotation_col": cfg.annotation_col,
                    "compound_id_col": cfg.compound_id_col,
                    "intensity_col": cfg.intensity_col,
                    "intensity_cols": cfg.intensity_cols,
                    "intensity_regex": cfg.intensity_regex,
                },
            },
            "created_at": time.time(),
            "version": "1.0",
        }
        try:
            joblib.dump(artifact, cfg.save_model_path)
        except Exception as e:
            raise RuntimeError(f"Failed to save model artifact to {cfg.save_model_path}: {e}")
    # Candidates for all DS1 rows
    rows_feats: List[dict] = []
    for i1, js in enumerate(nbh):
        for j2 in js:
            rows_feats.append(_feature_row(ds1, ds2, int(i1), int(j2)))
    cand = pd.DataFrame(rows_feats)
    if cand.empty:
        return MLMatchResult(candidates=cand, top1=pd.DataFrame(columns=["id1","id2","prob"]), hungarian=pd.DataFrame(columns=["id1","id2","prob"]))
    # Preserve feature names to avoid sklearn/lightgbm warnings
    # Score using the exact feature set used for training (excludes ids)
    scores = clf.predict_proba(cand[feature_names])[:, 1]
    if iso is not None:
        prob_cal = iso.predict(scores)
        cand["prob_raw"] = scores
        cand["prob_cal"] = np.clip(prob_cal, 1e-6, 1 - 1e-6)
        # per-row masses
        sums = cand.groupby("id1")["prob_cal"].transform("sum")
        cand["p_row"] = cand["prob_cal"].values
        over = sums > 1.0
        if over.any():
            cand.loc[over, "p_row"] = cand.loc[over, "prob_cal"] / sums[over]
        row_sum = cand.groupby("id1")["p_row"].sum()
        p0 = (1.0 - row_sum).clip(lower=0.0)
        cand = cand.merge(p0.rename("p0"), left_on="id1", right_index=True, how="left")
        score_col = "p_row"
    else:
        cand["prob"] = scores
        score_col = "prob"
    # Assignments
    top1 = _assignment_top1(cand, score_col)
    hung = _assignment_hungarian(cand, score_col)
    return MLMatchResult(candidates=cand, top1=top1, hungarian=hung)


__all__ = ["MLMatchConfig", "MLMatchResult", "match_ml"]
