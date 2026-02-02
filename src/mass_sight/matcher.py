"""massSight: probabilistic drift-aware optimal transport matcher.

This matcher builds drift-corrected residuals (RT, ppm, ROI, intensity),
scores candidates with a heavy-tailed likelihood, then applies entropic OT
to encourage global 1-1 consistency. It returns per-edge OT weights and
row-normalized probabilities for top-1 selection.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .structure_graph import (
    StructureGraph,
    build_structure_graph,
    build_structure_graph_local_rt,
    build_structure_graph_stable,
    coerce_expression,
)
from .lcms_utils import (
    ADDUCT_MASS_POS,
    ADDUCT_MASS_NEG,
    PROTON_MASS,
    SchemaConfig,
    infer_rt_unit_scale_to_minutes,
    normalize_schema,
)


@dataclass
class MassSightConfig:
    """Configuration for massSight probabilistic feature matching."""

    # Initial candidate windows
    ppm: float = 20.0
    # Cross-study default: disable hard RT candidate gating; handle RT via drift residuals.
    rt: float = 0.0
    # Intensity-based candidate gating is only applied when `use_intensity=True`.
    # Cross-study default: disable intensity-based filtering because intensities are not
    # comparable across instruments and preprocessing pipelines.
    intensity_log_diff_max: float = 5.0

    # Candidate expansion via discrete m/z shift hypotheses (experimental).
    #
    # We currently disable this by default because it degrades performance on the
    # canonical MW NEW_V1 benchmark. The functionality remains available for
    # targeted stress tests / future work.
    mz_shift_mode: str = "none"  # "none" | "auto" | "manual"
    # Positive-mode shift family selection (used only when polarity="positive").
    # - "conservative": ±H only
    # - "expanded": include common POS adducts (H/NH4/Na/K) in addition to ±H
    mz_shift_pos_adducts: str = "conservative"  # "conservative" | "expanded"
    # Negative-mode shift family selection (used only when polarity="negative").
    # - "conservative": ±H only
    # - "expanded": include common NEG adducts (Cl/formate/acetate) in addition to ±H
    mz_shift_neg_adducts: str = "conservative"  # "conservative" | "expanded"
    mz_shift_top_k: int = 1
    mz_shift_min_support_count: int = 10
    mz_shift_min_support_frac: float = 0.01
    mz_shift_max_abs_da: float = 40.0
    mz_shift_manual_deltas_da: List[float] = field(default_factory=list)
    mz_shift_include_isotope: bool = True
    mz_shift_penalty: float = 0.5  # set to 0 to disable
    # Auto-shift selection: require a shift to be a strong outlier in support counts.
    # 0 disables this filter (legacy behavior).
    mz_shift_min_robust_z: float = 0.0
    # Per-feature gate: only allow nonzero shifts for study_a features with no delta=0
    # candidate within the tight_ppm window (reduces spurious expansion on dense datasets).
    mz_shift_per_feature_tight_ppm_gate: bool = True

    # Candidate RT gating: optionally enable gating but disable it automatically when
    # the estimated RT drift (from m/z-tight anchors) is too large.
    rt_gate_mode: str = "none"  # "none" | "adaptive-disable"
    rt_gate_disable_p95: float = 1.0
    rt_gate_anchor_ppm: Optional[float] = None
    rt_gate_min_count: int = 25

    # Retention handling under missing/degenerate RT.
    # Some untargeted datasets (e.g., certain MW `untarg_data` exports) have no meaningful RT
    # in their feature IDs. In these cases RT/ROI residuals can become non-finite and collapse
    # OT to the null. In `auto` mode, massSight disables retention terms when RT quality is poor.
    retention_mode: str = "auto"  # "auto" | "use" | "mz_only"
    retention_min_finite_frac: float = 0.5
    retention_min_unique: int = 20
    retention_min_unique_frac: float = 0.01
    retention_round_decimals: int = 3

    # RT unit normalization + per-pair RT transferability.
    #
    # Some public feature tables store RT in seconds while others store minutes. In addition,
    # even minutes-scale RT is not always transferable across studies (e.g., different HILIC
    # columns/gradients). In `auto` mode we use m/z-tight mutual-NN anchors to decide whether
    # to (a) fit a simple linear RT transfer (scale+offset) or (b) ignore RT/ROI terms.
    normalize_rt_units: bool = True
    rt_transferability: str = "auto"  # "off" | "auto"
    rt_transfer_anchor_ppm: float = 5.0
    rt_transfer_min_anchors: int = 25
    rt_transfer_spearman_min: float = 0.8
    rt_transfer_p95_abs_resid_max: float = 2.0
    rt_transfer_max_anchors: int = 5000
    rt_transfer_ransac_iters: int = 2000

    # Tight windows for drift estimation
    tight_ppm: float = 7.0
    tight_rt: float = 0.0
    # Optional retention-order gate for tight anchors (in [0,1] since `ro` is a percentile rank).
    # This can be more stable than absolute RT windows across labs/instruments.
    tight_roi: float = 0.0
    # Strategy for selecting tight anchors used for drift fitting:
    #  - "all": keep all within tight windows (legacy)
    #  - "nn": per-feature nearest neighbor (by |ppm|) within tight windows
    #  - "mnn": mutual nearest neighbors (higher precision, fewer anchors)
    tight_anchor_strategy: str = "mnn"
    # Intensity-based tight-anchor filtering is only applied when `use_intensity=True`.
    tight_intensity_log_diff: float = 5.0

    # Spline smoothing parameters
    rt_spline_df: int = 10
    ppm_spline_df: int = 10
    roi_spline_df: int = 10
    rt_drift_x: str = "rt"  # "rt" | "ro" (retention order)
    # Cross-study default: disable explicit RT drift correction (often overfits sparse anchors).
    rt_drift_model: str = "none"  # "none" | "lowess" | "linear" | "robust-linear"
    ppm_drift_model: str = "linear"  # "none" | "lowess" | "linear" | "robust-linear"
    # Retention-order drift correction model (delta_roi as a function of roi_x).
    # "none" disables explicit ROI drift correction (roi_fit=0 -> roi_resid=delta_roi).
    roi_drift_model: str = "none"  # "none" | "lowess" | "linear" | "robust-linear"

    # Drift uncertainty propagation (optional)
    # (1) Anchor-bootstrap drift: fit the *main* drift models on a bootstrap resample of tight matches.
    drift_fit_bootstrap: bool = False
    drift_bootstrap_seed: int = 0

    # Likelihood weights
    use_intensity: bool = False
    w_ppm: float = 1.0
    w_rt: float = 1.0
    w_roi: float = 1.0
    w_int: float = 0.0

    # Student-t likelihood parameters
    t_df: float = 3.0
    scale_floor: float = 1e-3
    use_anchor_scales: bool = True
    anchor_scale_min_count: int = 30

    # OT parameters
    # OT mode:
    # - "balanced": (default) entropic OT with row + column marginals (and optional null column)
    # - "unbalanced": entropic OT with relaxed KL-penalized marginals (better supports many-to-one)
    # - "semi_relaxed": constrain rows only (columns have infinite capacity; effectively a per-row softmax at temperature `ot_epsilon`)
    ot_mode: str = "balanced"  # "balanced" | "unbalanced" | "semi_relaxed"
    ot_epsilon: float = 0.5
    ot_max_iter: int = 200
    ot_tol: float = 1e-4
    ot_min_mass: float = 1e-12
    # Unbalanced OT (KL-relaxed marginals) parameters. Larger values enforce marginals more strongly.
    # Effective exponent: tau = lambda / (lambda + ot_epsilon), so lambda >> epsilon approximates balanced OT.
    ot_unbalanced_lambda_row: float = 10.0
    ot_unbalanced_lambda_col: float = 1.0
    # Cross-study default: permit explicit no-match via a null OT column.
    allow_unmatched: bool = True
    null_mass: float = 0.1
    null_loglik_quantile: float = 0.2
    null_loglik: Optional[float] = None
    max_p0: Optional[float] = None

    # Structure-aware matching (optional)
    use_structure: bool = False
    structure_alpha: float = 0.0
    structure_k: int = 20
    structure_corr: str = "spearman"
    structure_corr_estimator: str = "sample"  # "sample" | "ledoit_wolf"
    structure_abs: bool = True
    structure_min_samples: int = 25
    # Graph construction strategy
    structure_graph: str = "knn"  # "knn" | "stable" | "local_rt"
    # RT-local correlation graph parameters (used only when structure_graph="local_rt").
    structure_rt_window: float = 0.05
    structure_max_candidates: int = 200
    structure_n_boot: int = 10
    structure_subsample_frac: float = 0.8
    structure_min_freq: float = 0.6
    # Adaptive weighting: apply structure more strongly on ambiguous rows.
    structure_adaptive_alpha: bool = False
    structure_margin_tau: float = 0.2  # 0 disables adaptive scaling
    structure_margin_gamma: float = 1.0
    # Annealing across structure iterations.
    structure_anneal: bool = False
    structure_alpha_ramp: bool = True
    structure_ot_epsilon_mult_start: float = 1.0  # >1 starts with higher entropy
    structure_iters: int = 2
    structure_eps: float = 1e-8
    structure_bonus_mode: str = "mass"
    # Optional: when RT is deemed non-transferable (rt_policy="ignore") and expression matrices
    # are provided, automatically enable structure refinement with this alpha (>0).
    structure_alpha_on_nontransferable_rt: float = 0.0

    # Schema columns
    mz_col: str = "MZ"
    rt_col: str = "RT"
    annotation_col: str = "Annotation_ID"
    compound_id_col: str = "Compound_ID"
    intensity_col: Optional[str] = "Intensity"

    mz_col_study_a: Optional[str] = None
    mz_col_study_b: Optional[str] = None
    rt_col_study_a: Optional[str] = None
    rt_col_study_b: Optional[str] = None
    annotation_col_study_a: Optional[str] = None
    annotation_col_study_b: Optional[str] = None
    compound_id_col_study_a: Optional[str] = None
    compound_id_col_study_b: Optional[str] = None

    polarity: str = "positive"


@dataclass
class MatchResult:
    """Result of massSight matching."""

    candidates: pd.DataFrame
    top1: pd.DataFrame
    ot_summary: Dict[str, float]
    drift_params: Dict[str, object]


_C13_SHIFT_DA = 1.00335483507


def _compute_rt_quality(rt: np.ndarray, cfg: MassSightConfig) -> Dict[str, object]:
    rt = np.asarray(rt, dtype=float)
    n = int(rt.size)
    if n == 0:
        return {
            "rt_finite_frac": 0.0,
            "rt_n_finite": 0,
            "rt_unique": 0,
            "rt_unique_frac": 0.0,
            "rt_ok": False,
        }

    finite = np.isfinite(rt)
    n_finite = int(np.sum(finite))
    finite_frac = float(n_finite / n) if n else 0.0

    unique = 0
    unique_frac = 0.0
    if n_finite:
        round_decimals = int(getattr(cfg, "retention_round_decimals", 3) or 3)
        rt_round = np.round(rt[finite], round_decimals)
        unique = int(len(np.unique(rt_round)))
        unique_frac = float(unique / max(n_finite, 1))

    min_finite_frac = float(getattr(cfg, "retention_min_finite_frac", 0.5) or 0.5)
    min_unique = int(getattr(cfg, "retention_min_unique", 20) or 20)
    min_unique_frac = float(getattr(cfg, "retention_min_unique_frac", 0.01) or 0.01)

    ok = finite_frac >= min_finite_frac and (unique >= min_unique or unique_frac >= min_unique_frac)
    return {
        "rt_finite_frac": float(finite_frac),
        "rt_n_finite": int(n_finite),
        "rt_unique": int(unique),
        "rt_unique_frac": float(unique_frac),
        "rt_ok": bool(ok),
    }


def _maybe_disable_retention(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Tuple[MassSightConfig, Dict[str, object]]:
    mode_req = str(getattr(cfg, "retention_mode", "auto") or "auto").lower().strip()
    if mode_req in {"mz-only", "mz", "mz_only"}:
        mode_req = "mz_only"
    if mode_req not in {"auto", "use", "mz_only"}:
        raise ValueError(
            f"Unsupported retention_mode: {getattr(cfg, 'retention_mode', None)!r} "
            "(expected 'auto', 'use', or 'mz_only')."
        )

    rt_a = study_a["RT"].to_numpy(dtype=float) if "RT" in study_a.columns else np.array([], dtype=float)
    rt_b = study_b["RT"].to_numpy(dtype=float) if "RT" in study_b.columns else np.array([], dtype=float)
    q_a = _compute_rt_quality(rt_a, cfg)
    q_b = _compute_rt_quality(rt_b, cfg)
    rt_ok_a = bool(q_a.get("rt_ok", False))
    rt_ok_b = bool(q_b.get("rt_ok", False))

    disable = False
    reason: Optional[str] = None
    if mode_req == "mz_only":
        disable = True
        reason = "forced_mz_only"
    elif mode_req == "auto":
        if not rt_ok_a or not rt_ok_b:
            disable = True
            if not rt_ok_a and not rt_ok_b:
                reason = "rt_bad_both"
            elif not rt_ok_a:
                reason = "rt_bad_study_a"
            else:
                reason = "rt_bad_study_b"

    cfg_eff = cfg
    mode_eff = mode_req
    if disable:
        cfg_eff = replace(
            cfg_eff,
            w_rt=0.0,
            w_roi=0.0,
            rt=0.0,
            tight_rt=0.0,
            tight_roi=0.0,
            rt_drift_model="none",
            roi_drift_model="none",
        )
        mode_eff = "mz_only"

    diag: Dict[str, object] = {
        "retention_mode_requested": str(mode_req),
        "retention_mode_effective": str(mode_eff),
        "retention_disable_reason": str(reason) if reason else None,
        "rt_quality_study_a": q_a,
        "rt_quality_study_b": q_b,
        "rt_ok_study_a": bool(rt_ok_a),
        "rt_ok_study_b": bool(rt_ok_b),
    }
    return cfg_eff, diag


def _dedupe_deltas(deltas: List[float], *, tol: float = 1e-9) -> List[float]:
    out: List[float] = []
    for d in sorted([float(x) for x in deltas if np.isfinite(x)]):
        if out and abs(d - out[-1]) <= tol:
            continue
        out.append(d)
    return out


def _mz_shift_family(cfg: MassSightConfig) -> List[float]:
    polarity = str(cfg.polarity).lower().strip()
    if polarity == "positive":
        pos_mode = str(getattr(cfg, "mz_shift_pos_adducts", "conservative") or "conservative").lower().strip()
        if pos_mode not in {"conservative", "expanded"}:
            raise ValueError(f"Unsupported mz_shift_pos_adducts '{pos_mode}'.")
        if pos_mode == "conservative":
            adduct_masses = [0.0, float(PROTON_MASS)]
        else:
            adduct_masses = [0.0] + [float(x) for x in ADDUCT_MASS_POS.values()]
    elif polarity == "negative":
        neg_mode = str(getattr(cfg, "mz_shift_neg_adducts", "conservative") or "conservative").lower().strip()
        if neg_mode not in {"conservative", "expanded"}:
            raise ValueError(f"Unsupported mz_shift_neg_adducts '{neg_mode}'.")
        # Conservative NEG family: allow ±H to handle neutral vs [M-H]- mismatches.
        adduct_masses = [0.0, -float(PROTON_MASS)]
        if neg_mode == "expanded":
            # Expanded NEG family: include common adducts that drive large, discrete offsets.
            # This enables matching across studies that report different NEG adduct forms:
            #   [M-H]- ↔ [M+Cl]-   (Δ ≈ Cl+H)
            #   [M-H]- ↔ [M+FA-H]- (Δ ≈ formic acid)
            #   [M-H]- ↔ [M+Ac-H]- (Δ ≈ acetic acid)
            adduct_masses.extend(float(x) for x in ADDUCT_MASS_NEG.values())
    else:
        raise ValueError(f"Unsupported polarity '{cfg.polarity}'.")

    deltas = [a2 - a1 for a1 in adduct_masses for a2 in adduct_masses]
    if bool(cfg.mz_shift_include_isotope):
        deltas.extend([_C13_SHIFT_DA, -_C13_SHIFT_DA])

    max_abs = float(cfg.mz_shift_max_abs_da)
    # Expanded NEG deltas can exceed 40 Da (e.g., acetate), so ensure a usable cap.
    if polarity == "negative":
        neg_mode = str(getattr(cfg, "mz_shift_neg_adducts", "conservative") or "conservative").lower().strip()
        if neg_mode == "expanded" and (not np.isfinite(max_abs) or max_abs < 65.0):
            max_abs = 70.0
    if not np.isfinite(max_abs) or max_abs <= 0:
        max_abs = 40.0
    deltas = [d for d in deltas if np.isfinite(d) and abs(d) <= max_abs]

    deduped = _dedupe_deltas(deltas, tol=1e-9)
    if not deduped or abs(deduped[0]) > 1e-12:
        deduped = _dedupe_deltas([0.0] + deduped, tol=1e-9)
    return deduped


def _mz_shift_support_counts(
    mz1: np.ndarray,
    mz2_sorted: np.ndarray,
    *,
    ppm: float,
    deltas: List[float],
) -> Dict[float, int]:
    mz1 = np.asarray(mz1, dtype=float)
    mz2_sorted = np.asarray(mz2_sorted, dtype=float)

    ppm = float(ppm)
    if not np.isfinite(ppm) or ppm <= 0:
        ppm = 20.0

    counts: Dict[float, int] = {}
    for delta in deltas:
        d = float(delta)
        if not np.isfinite(d):
            continue
        support = 0
        for mz in mz1:
            target = float(mz) + d
            if not np.isfinite(target) or target <= 0:
                continue
            span = target * ppm / 1e6
            left = int(np.searchsorted(mz2_sorted, target - span, side="left"))
            right = int(np.searchsorted(mz2_sorted, target + span, side="right"))
            if left < right:
                support += 1
        counts[d] = int(support)
    return counts


def _select_mz_shifts(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Tuple[List[float], Dict[str, object]]:
    mode = str(getattr(cfg, "mz_shift_mode", "none") or "none").lower().strip()
    if mode not in {"none", "auto", "manual"}:
        raise ValueError(f"Unsupported mz_shift_mode: {cfg.mz_shift_mode!r}")

    diag: Dict[str, object] = {"mz_shift_mode": mode}
    if mode == "none":
        return [0.0], diag

    mz_b = study_b["MZ"].to_numpy(dtype=float)
    mz_b_sorted = np.sort(mz_b[np.isfinite(mz_b)])
    if mz_b_sorted.size == 0:
        return [0.0], diag

    if mode == "manual":
        manual = [float(x) for x in (cfg.mz_shift_manual_deltas_da or []) if np.isfinite(x)]
        deltas = _dedupe_deltas([0.0] + manual, tol=1e-9)
        diag["mz_shift_selected_deltas_da"] = deltas
        return deltas, diag

    family = _mz_shift_family(cfg)
    counts = _mz_shift_support_counts(
        study_a["MZ"].to_numpy(dtype=float),
        mz_b_sorted,
        ppm=float(cfg.ppm),
        deltas=family,
    )

    n1 = int(len(study_a))
    min_count = int(cfg.mz_shift_min_support_count)
    min_frac = float(cfg.mz_shift_min_support_frac)
    top_k = int(cfg.mz_shift_top_k)
    if top_k < 0:
        top_k = 0
    if not np.isfinite(min_frac) or min_frac < 0:
        min_frac = 0.0

    min_robust_z = float(getattr(cfg, "mz_shift_min_robust_z", 0.0) or 0.0)
    support_med = support_mad = float("nan")
    if min_robust_z > 0:
        nonzero_counts = np.array([int(v) for d, v in counts.items() if abs(float(d)) > 1e-12], dtype=float)
        if nonzero_counts.size > 0:
            support_med = float(np.median(nonzero_counts))
            support_mad = float(np.median(np.abs(nonzero_counts - support_med)))
            if not np.isfinite(support_mad) or support_mad <= 1e-12:
                support_mad = 1.0
    diag["mz_shift_min_robust_z"] = float(min_robust_z)
    diag["mz_shift_support_median_nonzero"] = float(support_med) if np.isfinite(support_med) else None
    diag["mz_shift_support_mad_nonzero"] = float(support_mad) if np.isfinite(support_mad) else None

    nonzero = []
    for d, c in counts.items():
        if abs(float(d)) <= 1e-12:
            continue
        frac = float(c) / float(n1) if n1 else 0.0
        if c >= min_count and frac >= min_frac:
            if min_robust_z > 0 and np.isfinite(support_med) and np.isfinite(support_mad):
                z = (float(c) - support_med) / float(support_mad)
                if z < min_robust_z:
                    continue
            nonzero.append((int(c), float(d)))
    nonzero.sort(key=lambda x: (-x[0], abs(x[1]), x[1]))

    selected = [0.0] + [d for _, d in nonzero[:top_k]]

    diag["mz_shift_selected_deltas_da"] = selected
    diag["mz_shift_top_k"] = int(cfg.mz_shift_top_k)
    diag["mz_shift_min_support_count"] = int(cfg.mz_shift_min_support_count)
    diag["mz_shift_min_support_frac"] = float(cfg.mz_shift_min_support_frac)
    diag["mz_shift_support_counts"] = {f"{float(d):+.12f}": int(c) for d, c in counts.items()}
    return selected, diag


def _candidate_neighborhood_indices(
    mz2_sorted: np.ndarray,
    order: np.ndarray,
    *,
    mz_target: float,
    ppm: float,
) -> np.ndarray:
    span = mz_target * ppm / 1e6
    lo, hi = mz_target - span, mz_target + span
    left = np.searchsorted(mz2_sorted, lo, side="left")
    right = np.searchsorted(mz2_sorted, hi, side="right")
    return order[left:right]


def _ppm_signed_over_mean(mz_x: float, mz_y: float) -> float:
    denom = max(0.5 * (float(mz_x) + float(mz_y)), 1e-12)
    return (float(mz_y) - float(mz_x)) / denom * 1e6


def _nearest_neighbor_mz_matches(
    mz_query: np.ndarray,
    mz_target: np.ndarray,
    *,
    mz_ppm: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    For each mz_query[i], find the nearest mz_target[j] within `mz_ppm` (ppm over mean m/z).

    Returns:
      - best_idx: shape (n_query,), -1 when no match within tolerance
      - delta_ppm: signed ppm error for the chosen neighbor, NaN when no match
    """
    mz_query = np.asarray(mz_query, dtype=float)
    mz_target = np.asarray(mz_target, dtype=float)

    n_q = int(mz_query.size)
    best_idx = np.full(n_q, -1, dtype=int)
    delta_ppm = np.full(n_q, np.nan, dtype=float)
    if n_q == 0 or int(mz_target.size) == 0:
        return best_idx, delta_ppm

    mask_t = np.isfinite(mz_target) & (mz_target > 0)
    if int(mask_t.sum()) == 0:
        return best_idx, delta_ppm

    order = np.argsort(mz_target[mask_t])
    idx_t = np.nonzero(mask_t)[0][order]
    mz_sorted = mz_target[idx_t]
    n_t = int(mz_sorted.size)
    if n_t == 0:
        return best_idx, delta_ppm

    mask_q = np.isfinite(mz_query) & (mz_query > 0)
    if int(mask_q.sum()) == 0:
        return best_idx, delta_ppm

    q = mz_query[mask_q]
    pos = np.searchsorted(mz_sorted, q)
    pos0 = np.clip(pos, 0, n_t - 1)
    pos1 = np.clip(pos - 1, 0, n_t - 1)

    cand0 = idx_t[pos0]
    cand1 = idx_t[pos1]
    diff0 = np.abs(mz_target[cand0] - q)
    diff1 = np.abs(mz_target[cand1] - q)
    use0 = diff0 <= diff1
    cand = np.where(use0, cand0, cand1).astype(int)

    denom = 0.5 * (mz_target[cand] + q)
    denom = np.where(denom <= 1e-12, 1e-12, denom)
    dppm = (mz_target[cand] - q) / denom * 1e6
    ok = np.abs(dppm) <= float(mz_ppm)

    q_idx = np.nonzero(mask_q)[0]
    best_idx[q_idx[ok]] = cand[ok]
    delta_ppm[q_idx[ok]] = dppm[ok]
    return best_idx, delta_ppm


def _infer_pair_rt_transferability(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Dict[str, object]:
    """
    Decide whether RT is transferable across a dataset pair and, if so, fit a linear RT transform.

    Uses m/z-tight mutual-nearest-neighbor (MNN) anchors as a high-precision proxy for shared ions.

    Returns a dict with:
      - rt_policy: 'linear' | 'ignore' | 'unknown'
      - rt_scale, rt_offset: mapping from dataset B -> dataset A (rt_a ≈ scale*rt_b + offset)
      - anchor diagnostics: n_anchors, spearman_rt, p95_abs_resid_rt, inlier_frac, etc.
    """
    mz1 = pd.to_numeric(study_a.get("MZ"), errors="coerce").to_numpy(dtype=float, copy=False)
    mz2 = pd.to_numeric(study_b.get("MZ"), errors="coerce").to_numpy(dtype=float, copy=False)
    rt1 = pd.to_numeric(study_a.get("RT"), errors="coerce").to_numpy(dtype=float, copy=False)
    rt2 = pd.to_numeric(study_b.get("RT"), errors="coerce").to_numpy(dtype=float, copy=False)

    anchor_ppm = float(getattr(cfg, "rt_transfer_anchor_ppm", 5.0) or 5.0)
    min_anchors = int(getattr(cfg, "rt_transfer_min_anchors", 25) or 25)
    spearman_min = float(getattr(cfg, "rt_transfer_spearman_min", 0.8) or 0.8)
    p95_abs_resid_max = float(getattr(cfg, "rt_transfer_p95_abs_resid_max", 2.0) or 2.0)
    max_anchors = int(getattr(cfg, "rt_transfer_max_anchors", 5000) or 5000)
    max_anchors = max(max_anchors, 0)
    ransac_cap = int(getattr(cfg, "rt_transfer_ransac_iters", 2000) or 2000)

    out: Dict[str, object] = {
        "rt_policy": "unknown",
        "rt_scale": 1.0,
        "rt_offset": 0.0,
        "rt_anchor_ppm": float(anchor_ppm),
        "rt_min_anchors": int(min_anchors),
        "rt_n_anchors": 0,
        "rt_n_inliers": 0,
        "rt_inlier_frac": float("nan"),
        "rt_ransac_tol": float("nan"),
        "rt_spearman": float("nan"),
        "rt_median_abs_resid": float("nan"),
        "rt_p95_abs_resid": float("nan"),
        "rt_note": "",
    }

    if int(mz1.size) == 0 or int(mz2.size) == 0:
        out["rt_note"] = "empty_mz"
        return out

    best2_for_1, _ = _nearest_neighbor_mz_matches(mz1, mz2, mz_ppm=float(anchor_ppm))
    best1_for_2, _ = _nearest_neighbor_mz_matches(mz2, mz1, mz_ppm=float(anchor_ppm))

    i = np.arange(mz1.size, dtype=int)
    j = best2_for_1
    mnn = (j >= 0) & (best1_for_2[j] == i)
    if not np.any(mnn):
        out["rt_note"] = "no_mnn_anchors"
        return out

    i = i[mnn]
    j = j[mnn]

    rt1_a = rt1[i]
    rt2_a = rt2[j]
    keep = np.isfinite(rt1_a) & np.isfinite(rt2_a)
    rt1_a = rt1_a[keep]
    rt2_a = rt2_a[keep]

    n = int(rt1_a.size)
    out["rt_n_anchors"] = int(n)
    if n < int(min_anchors):
        out["rt_note"] = "too_few_anchors"
        return out

    if max_anchors > 0 and n > max_anchors:
        rng = np.random.default_rng(0)
        sel = rng.choice(n, size=max_anchors, replace=False)
        rt1_a = rt1_a[sel]
        rt2_a = rt2_a[sel]
        n = int(rt1_a.size)
        out["rt_n_anchors"] = int(n)
        out["rt_note"] = "subsampled_anchors"

    # Choose an inlier tolerance on a minutes-like scale; proportional to RT range but bounded.
    with np.errstate(invalid="ignore"):
        r1_range = float(np.nanmax(rt1_a) - np.nanmin(rt1_a))
        r2_range = float(np.nanmax(rt2_a) - np.nanmin(rt2_a))
    r_range = float(min(r1_range, r2_range))
    tol = float(min(1.0, max(0.2, 0.05 * r_range)))
    out["rt_ransac_tol"] = float(tol)

    rng = np.random.default_rng(0)
    iters = int(min(max(ransac_cap, 0), max(200, 2 * n)))
    idx = np.arange(n, dtype=int)

    best_inliers = None
    best_count = -1
    best_scale = 1.0
    best_offset = 0.0
    for _ in range(iters):
        a, b = rng.choice(idx, size=2, replace=False)
        x1 = float(rt2_a[a])
        x2 = float(rt2_a[b])
        if not np.isfinite(x1) or not np.isfinite(x2) or abs(x2 - x1) <= 1e-6:
            continue
        y1 = float(rt1_a[a])
        y2 = float(rt1_a[b])
        if not np.isfinite(y1) or not np.isfinite(y2):
            continue
        scale = (y2 - y1) / (x2 - x1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        offset = y1 - scale * x1
        resid = (rt2_a * float(scale) + float(offset)) - rt1_a
        inliers = np.abs(resid) <= tol
        count = int(np.sum(inliers))
        if count > best_count:
            best_count = count
            best_scale = float(scale)
            best_offset = float(offset)
            best_inliers = inliers

    if best_inliers is None or best_count <= 0:
        out["rt_note"] = "ransac_failed"
        return out

    # Refine using least squares on inliers.
    in_mask = best_inliers
    rt1_in = rt1_a[in_mask]
    rt2_in = rt2_a[in_mask]
    out["rt_n_inliers"] = int(rt1_in.size)
    out["rt_inlier_frac"] = float(rt1_in.size / max(n, 1))

    X = np.vstack([rt2_in, np.ones_like(rt2_in)]).T
    beta = np.linalg.lstsq(X, rt1_in, rcond=None)[0]
    scale = float(beta[0])
    offset = float(beta[1])
    if not np.isfinite(scale) or scale <= 0:
        out["rt_note"] = "invalid_scale_refit"
        return out
    scale = float(np.clip(scale, 1.0 / 120.0, 120.0))

    rt2_hat_in = rt2_in * scale + offset
    resid_in = rt2_hat_in - rt1_in
    abs_resid_in = np.abs(resid_in)
    out["rt_scale"] = float(scale)
    out["rt_offset"] = float(offset)
    out["rt_median_abs_resid"] = float(np.median(abs_resid_in))
    out["rt_p95_abs_resid"] = float(np.quantile(abs_resid_in, 0.95))

    if rt1_in.size >= 3:
        r1 = pd.Series(rt1_in).rank(method="average").to_numpy(dtype=float)
        r2 = pd.Series(rt2_hat_in).rank(method="average").to_numpy(dtype=float)
        if float(np.std(r1)) > 0 and float(np.std(r2)) > 0:
            out["rt_spearman"] = float(np.corrcoef(r1, r2)[0, 1])

    # Decision rule: trust RT when there is a substantial inlier population supporting a monotone linear map.
    min_inlier_frac = 0.1
    transferable = (
        (out["rt_n_inliers"] >= int(min_anchors))
        and np.isfinite(out["rt_spearman"])
        and (float(out["rt_spearman"]) >= float(spearman_min))
        and (float(out["rt_p95_abs_resid"]) <= float(p95_abs_resid_max))
        and (float(out["rt_inlier_frac"]) >= float(min_inlier_frac))
    )
    if transferable:
        out["rt_policy"] = "linear"
        out["rt_note"] = "transferable"
    else:
        out["rt_policy"] = "ignore"
        out["rt_note"] = "nontransferable"
    return out


def _build_candidate_pairs_shifted(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
    shifts_da: List[float],
) -> pd.DataFrame:
    mz2 = study_b["MZ"].to_numpy(dtype=float)
    rt2 = study_b["RT"].to_numpy(dtype=float)
    order = np.argsort(mz2)
    mz2_sorted = mz2[order]
    rt2_sorted = rt2[order]
    mz1_all = study_a["MZ"].to_numpy(dtype=float)
    rt1_all = study_a["RT"].to_numpy(dtype=float)

    rt_gate = float(cfg.rt) if cfg.rt is not None else 0.0
    ppm = float(cfg.ppm)
    per_feature_gate = bool(getattr(cfg, "mz_shift_per_feature_tight_ppm_gate", False))
    tight_ppm_gate = float(cfg.tight_ppm)

    rows: List[dict] = []
    use_intensity_gate = bool(getattr(cfg, "use_intensity", False))
    for i1, (mz1, rt1) in enumerate(zip(mz1_all, rt1_all)):
        if not np.isfinite(mz1) or mz1 <= 0:
            continue
        shifts_iter = shifts_da
        if per_feature_gate and len(shifts_da) > 1 and np.isfinite(tight_ppm_gate) and tight_ppm_gate > 0:
            js0 = _candidate_neighborhood_indices(mz2_sorted, order, mz_target=float(mz1), ppm=tight_ppm_gate)
            if js0.size > 0:
                shifts_iter = [0.0]
        best_by_j: Dict[int, Tuple[float, float]] = {}
        for delta in shifts_iter:
            d = float(delta)
            mz_target = float(mz1) + d
            if not np.isfinite(mz_target) or mz_target <= 0:
                continue
            js = _candidate_neighborhood_indices(mz2_sorted, order, mz_target=mz_target, ppm=ppm)
            if js.size == 0:
                continue
            for j2 in js.tolist():
                if rt_gate > 0:
                    rt2_val = float(study_b.at[int(j2), "RT"])
                    if abs(rt2_val - float(rt1)) > rt_gate:
                        continue
                mz2_val = float(study_b.at[int(j2), "MZ"])
                abs_ppm = abs(_ppm_signed_over_mean(mz_target, mz2_val))
                prev = best_by_j.get(int(j2))
                if prev is None or abs_ppm < prev[0] - 1e-12 or (
                    abs(abs_ppm - prev[0]) <= 1e-12 and (abs(d), d) < (abs(prev[1]), prev[1])
                ):
                    best_by_j[int(j2)] = (float(abs_ppm), d)

        if not best_by_j:
            continue

        for j2, (_, d_best) in best_by_j.items():
            mz2_val = float(study_b.at[int(j2), "MZ"])
            rt2_val = float(study_b.at[int(j2), "RT"])
            ro1 = float(study_a.at[int(i1), "ro"]) if "ro" in study_a.columns else 0.0
            ro2 = float(study_b.at[int(j2), "ro"]) if "ro" in study_b.columns else 0.0
            il1 = float(study_a.at[int(i1), "Intensity_log10"]) if "Intensity_log10" in study_a.columns else 0.0
            il2 = float(study_b.at[int(j2), "Intensity_log10"]) if "Intensity_log10" in study_b.columns else 0.0
            feat = {
                "id1": int(i1),
                "id2": int(j2),
                "mz1": float(mz1),
                "mz2": float(mz2_val),
                "rt1": float(rt1),
                "rt2": float(rt2_val),
                "ro1": float(ro1),
                "ro2": float(ro2),
                "int1": float(il1),
                "int2": float(il2),
            }
            feat["mz_shift_da"] = float(d_best)
            feat["intensity_log_diff"] = il2 - il1

            feat["roi_x"] = float(ro1)
            feat["roi_y"] = float(ro2)
            feat["delta_roi"] = feat["roi_y"] - feat["roi_x"]

            if use_intensity_gate:
                if abs(feat["intensity_log_diff"]) <= float(cfg.intensity_log_diff_max):
                    rows.append(feat)
            else:
                rows.append(feat)

    return pd.DataFrame(rows)


def _build_tight_matches_shifted(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
    shifts_da: List[float],
) -> pd.DataFrame:
    study_a = study_a.reset_index(drop=True)
    study_b = study_b.reset_index(drop=True)

    mz2 = study_b["MZ"].to_numpy(dtype=float)
    rt2 = study_b["RT"].to_numpy(dtype=float)
    ro2 = study_b["ro"].to_numpy(dtype=float) if "ro" in study_b.columns else np.zeros(len(study_b), dtype=float)
    il2 = (
        study_b["Intensity_log10"].to_numpy(dtype=float)
        if "Intensity_log10" in study_b.columns
        else np.zeros(len(study_b), dtype=float)
    )

    order2 = np.argsort(mz2)
    mz2_sorted = mz2[order2]

    mz1 = study_a["MZ"].to_numpy(dtype=float)
    rt1 = study_a["RT"].to_numpy(dtype=float)
    ro1 = study_a["ro"].to_numpy(dtype=float) if "ro" in study_a.columns else np.zeros(len(study_a), dtype=float)
    il1 = (
        study_a["Intensity_log10"].to_numpy(dtype=float)
        if "Intensity_log10" in study_a.columns
        else np.zeros(len(study_a), dtype=float)
    )

    tight_ppm = float(getattr(cfg, "tight_ppm", 0.0) or 0.0)
    tight_rt = float(getattr(cfg, "tight_rt", 0.0) or 0.0)
    tight_roi = float(getattr(cfg, "tight_roi", 0.0) or 0.0)
    use_intensity_gate = bool(getattr(cfg, "use_intensity", False))
    tight_int_gate = float(getattr(cfg, "tight_intensity_log_diff", 0.0) or 0.0)

    strategy = str(getattr(cfg, "tight_anchor_strategy", "all") or "all").lower().strip()
    if strategy in {"nearest", "nearest-neighbor", "nearest_neighbor"}:
        strategy = "nn"
    if strategy in {"mutual", "mutual-nearest", "mutual_nn", "mutual-nearest-neighbor"}:
        strategy = "mnn"
    if strategy not in {"all", "nn", "mnn"}:
        raise ValueError(
            f"Unsupported tight_anchor_strategy: {getattr(cfg, 'tight_anchor_strategy', None)!r} "
            "(expected 'all', 'nn', or 'mnn')."
        )

    if not np.isfinite(tight_ppm) or tight_ppm <= 0:
        return pd.DataFrame(
            columns=[
                "id1",
                "id2",
                "mz_x",
                "mz_y",
                "rt_x",
                "rt_y",
                "roi_x",
                "roi_y",
                "int_x",
                "int_y",
                "delta_ppm",
                "delta_rt",
                "delta_roi",
                "intensity_log_diff",
                "mz_shift_da",
            ]
        )

    shifts_da = [float(d) for d in shifts_da] if shifts_da else [0.0]
    per_feature_gate = bool(getattr(cfg, "mz_shift_per_feature_tight_ppm_gate", False))

    def _passes_gates(i: int, j: int) -> bool:
        if tight_rt > 0 and abs(float(rt2[j]) - float(rt1[i])) > tight_rt:
            return False
        if tight_roi > 0 and abs(float(ro2[j]) - float(ro1[i])) > tight_roi:
            return False
        if use_intensity_gate and abs(float(il2[j]) - float(il1[i])) > tight_int_gate:
            return False
        return True

    rows: List[dict] = []

    if strategy == "all":
        for i in range(len(study_a)):
            mz1_val = float(mz1[i])
            if not np.isfinite(mz1_val) or mz1_val <= 0:
                continue

            shifts_iter = shifts_da
            if per_feature_gate and len(shifts_da) > 1:
                span0 = mz1_val * tight_ppm / 1e6
                left0 = int(np.searchsorted(mz2_sorted, mz1_val - span0, side="left"))
                right0 = int(np.searchsorted(mz2_sorted, mz1_val + span0, side="right"))
                if left0 < right0:
                    shifts_iter = [0.0]

            best_by_j: Dict[int, Tuple[float, float]] = {}
            for d in shifts_iter:
                mz_target = mz1_val + float(d)
                if not np.isfinite(mz_target) or mz_target <= 0:
                    continue
                span = mz_target * tight_ppm / 1e6
                left = int(np.searchsorted(mz2_sorted, mz_target - span, side="left"))
                right = int(np.searchsorted(mz2_sorted, mz_target + span, side="right"))
                if left >= right:
                    continue
                for idx in range(left, right):
                    j = int(order2[idx])
                    if not _passes_gates(i, j):
                        continue
                    mz2_val = float(mz2[j])
                    abs_ppm = abs(_ppm_signed_over_mean(mz_target, mz2_val))
                    prev = best_by_j.get(j)
                    if prev is None or abs_ppm < prev[0] - 1e-12 or (
                        abs(abs_ppm - prev[0]) <= 1e-12 and (abs(float(d)), float(d)) < (abs(prev[1]), prev[1])
                    ):
                        best_by_j[j] = (float(abs_ppm), float(d))

            for j, (_, d_best) in best_by_j.items():
                mz_target = mz1_val + float(d_best)
                mz2_val = float(mz2[j])
                delta_ppm = _ppm_signed_over_mean(mz_target, mz2_val)
                delta_rt = float(rt2[j]) - float(rt1[i])
                delta_roi = float(ro2[j]) - float(ro1[i])
                int_diff = float(il2[j]) - float(il1[i])

                rows.append(
                    {
                        "id1": int(i),
                        "id2": int(j),
                        "mz_x": float(mz_target),
                        "mz_y": float(mz2_val),
                        "rt_x": float(rt1[i]),
                        "rt_y": float(rt2[j]),
                        "roi_x": float(ro1[i]),
                        "roi_y": float(ro2[j]),
                        "int_x": float(il1[i]),
                        "int_y": float(il2[j]),
                        "delta_ppm": float(delta_ppm),
                        "delta_rt": float(delta_rt),
                        "delta_roi": float(delta_roi),
                        "intensity_log_diff": float(int_diff),
                        "mz_shift_da": float(d_best),
                    }
                )

        return pd.DataFrame(rows)

    def _best_match_a_to_b() -> Tuple[np.ndarray, np.ndarray]:
        best_j = np.full(len(study_a), -1, dtype=int)
        best_shift = np.zeros(len(study_a), dtype=float)
        best_abs_ppm = np.full(len(study_a), np.inf, dtype=float)
        best_abs_shift = np.full(len(study_a), np.inf, dtype=float)
        best_abs_drt = np.full(len(study_a), np.inf, dtype=float)

        for i in range(len(study_a)):
            mz1_val = float(mz1[i])
            if not np.isfinite(mz1_val) or mz1_val <= 0:
                continue

            shifts_iter = shifts_da
            if per_feature_gate and len(shifts_da) > 1:
                span0 = mz1_val * tight_ppm / 1e6
                left0 = int(np.searchsorted(mz2_sorted, mz1_val - span0, side="left"))
                right0 = int(np.searchsorted(mz2_sorted, mz1_val + span0, side="right"))
                if left0 < right0:
                    shifts_iter = [0.0]

            for d in shifts_iter:
                mz_target = mz1_val + float(d)
                if not np.isfinite(mz_target) or mz_target <= 0:
                    continue
                span = mz_target * tight_ppm / 1e6
                left = int(np.searchsorted(mz2_sorted, mz_target - span, side="left"))
                right = int(np.searchsorted(mz2_sorted, mz_target + span, side="right"))
                if left >= right:
                    continue

                for idx in range(left, right):
                    j = int(order2[idx])
                    if not _passes_gates(i, j):
                        continue
                    mz2_val = float(mz2[j])
                    abs_ppm = abs(_ppm_signed_over_mean(mz_target, mz2_val))
                    if not np.isfinite(abs_ppm):
                        continue
                    abs_shift = abs(float(d))
                    abs_drt = abs(float(rt2[j]) - float(rt1[i]))
                    if (
                        abs_ppm < float(best_abs_ppm[i]) - 1e-12
                        or (
                            abs(abs_ppm - float(best_abs_ppm[i])) <= 1e-12
                            and (abs_shift, abs_drt) < (float(best_abs_shift[i]), float(best_abs_drt[i]))
                        )
                    ):
                        best_abs_ppm[i] = float(abs_ppm)
                        best_abs_shift[i] = float(abs_shift)
                        best_abs_drt[i] = float(abs_drt)
                        best_j[i] = int(j)
                        best_shift[i] = float(d)

        return best_j, best_shift

    if strategy == "nn":
        best_j, best_shift = _best_match_a_to_b()
        for i in range(len(study_a)):
            j = int(best_j[i])
            if j < 0:
                continue
            d_best = float(best_shift[i])
            mz_target = float(mz1[i]) + d_best
            mz2_val = float(mz2[j])
            delta_ppm = _ppm_signed_over_mean(mz_target, mz2_val)
            delta_rt = float(rt2[j]) - float(rt1[i])
            delta_roi = float(ro2[j]) - float(ro1[i])
            int_diff = float(il2[j]) - float(il1[i])
            rows.append(
                {
                    "id1": int(i),
                    "id2": int(j),
                    "mz_x": float(mz_target),
                    "mz_y": float(mz2_val),
                    "rt_x": float(rt1[i]),
                    "rt_y": float(rt2[j]),
                    "roi_x": float(ro1[i]),
                    "roi_y": float(ro2[j]),
                    "int_x": float(il1[i]),
                    "int_y": float(il2[j]),
                    "delta_ppm": float(delta_ppm),
                    "delta_rt": float(delta_rt),
                    "delta_roi": float(delta_roi),
                    "intensity_log_diff": float(int_diff),
                    "mz_shift_da": float(d_best),
                }
            )
        return pd.DataFrame(rows)

    # mutual nearest neighbors
    best_j, best_shift = _best_match_a_to_b()

    mz1_shift_sorted: List[np.ndarray] = []
    order1_by_shift: List[np.ndarray] = []
    for d in shifts_da:
        mz1s = mz1 + float(d)
        order1 = np.argsort(mz1s)
        mz1_shift_sorted.append(mz1s[order1])
        order1_by_shift.append(order1)

    best_i_rev = np.full(len(study_b), -1, dtype=int)
    best_shift_rev = np.zeros(len(study_b), dtype=float)
    best_abs_ppm_rev = np.full(len(study_b), np.inf, dtype=float)
    best_abs_shift_rev = np.full(len(study_b), np.inf, dtype=float)
    best_abs_drt_rev = np.full(len(study_b), np.inf, dtype=float)

    for j in range(len(study_b)):
        mz2_val = float(mz2[j])
        if not np.isfinite(mz2_val) or mz2_val <= 0:
            continue
        span2 = mz2_val * tight_ppm / 1e6
        for k, d in enumerate(shifts_da):
            mz1s_sorted = mz1_shift_sorted[k]
            order1 = order1_by_shift[k]
            left = int(np.searchsorted(mz1s_sorted, mz2_val - span2, side="left"))
            right = int(np.searchsorted(mz1s_sorted, mz2_val + span2, side="right"))
            if left >= right:
                continue
            for idx in range(left, right):
                i = int(order1[idx])
                if not _passes_gates(i, j):
                    continue
                mz_target = float(mz1[i]) + float(d)
                abs_ppm = abs(_ppm_signed_over_mean(mz_target, mz2_val))
                if not np.isfinite(abs_ppm):
                    continue
                abs_shift = abs(float(d))
                abs_drt = abs(float(rt2[j]) - float(rt1[i]))
                if (
                    abs_ppm < float(best_abs_ppm_rev[j]) - 1e-12
                    or (
                        abs(abs_ppm - float(best_abs_ppm_rev[j])) <= 1e-12
                        and (abs_shift, abs_drt) < (float(best_abs_shift_rev[j]), float(best_abs_drt_rev[j]))
                    )
                ):
                    best_abs_ppm_rev[j] = float(abs_ppm)
                    best_abs_shift_rev[j] = float(abs_shift)
                    best_abs_drt_rev[j] = float(abs_drt)
                    best_i_rev[j] = int(i)
                    best_shift_rev[j] = float(d)

    for i in range(len(study_a)):
        j = int(best_j[i])
        if j < 0:
            continue
        if int(best_i_rev[j]) != int(i):
            continue
        d_best = float(best_shift[i])
        mz_target = float(mz1[i]) + d_best
        mz2_val = float(mz2[j])
        delta_ppm = _ppm_signed_over_mean(mz_target, mz2_val)
        delta_rt = float(rt2[j]) - float(rt1[i])
        delta_roi = float(ro2[j]) - float(ro1[i])
        int_diff = float(il2[j]) - float(il1[i])
        rows.append(
            {
                "id1": int(i),
                "id2": int(j),
                "mz_x": float(mz_target),
                "mz_y": float(mz2_val),
                "rt_x": float(rt1[i]),
                "rt_y": float(rt2[j]),
                "roi_x": float(ro1[i]),
                "roi_y": float(ro2[j]),
                "int_x": float(il1[i]),
                "int_y": float(il2[j]),
                "delta_ppm": float(delta_ppm),
                "delta_rt": float(delta_rt),
                "delta_roi": float(delta_roi),
                "intensity_log_diff": float(int_diff),
                "mz_shift_da": float(d_best),
            }
        )

    return pd.DataFrame(rows)


def _estimate_rt_drift_p95_mz_tight(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    *,
    ppm: float,
    min_count: int,
) -> Tuple[float, int]:
    mz1 = study_a["MZ"].to_numpy(dtype=float)
    rt1 = study_a["RT"].to_numpy(dtype=float)
    mz2 = study_b["MZ"].to_numpy(dtype=float)
    rt2 = study_b["RT"].to_numpy(dtype=float)

    mask2 = np.isfinite(mz2)
    if int(mask2.sum()) == 0:
        return float("nan"), 0

    order = np.argsort(mz2[mask2])
    idx2 = np.nonzero(mask2)[0][order]
    mz2_sorted = mz2[idx2]

    deltas: List[float] = []
    for mz, rt in zip(mz1, rt1, strict=True):
        if not np.isfinite(mz) or mz <= 0:
            continue
        span = float(mz) * float(ppm) / 1e6
        left = int(np.searchsorted(mz2_sorted, float(mz) - span, side="left"))
        right = int(np.searchsorted(mz2_sorted, float(mz) + span, side="right"))
        if left >= right:
            continue
        cand = idx2[left:right]
        mz2_cand = mz2[cand]
        denom = 0.5 * (mz2_cand + float(mz))
        denom = np.where(denom <= 1e-12, 1e-12, denom)
        ppm_err = 1e6 * (mz2_cand - float(mz)) / denom
        j = int(cand[int(np.argmin(np.abs(ppm_err)))])
        deltas.append(abs(float(rt2[j]) - float(rt)))

    if len(deltas) < int(min_count):
        return float("nan"), int(len(deltas))
    return float(np.quantile(np.asarray(deltas, dtype=float), 0.95)), int(len(deltas))


def _fit_lowess_1d(x: np.ndarray, y: np.ndarray, frac: float = 0.2) -> np.ndarray:
    """Fit a LOWESS smoother (local linear regression with tricube weights)."""
    n = len(x)
    if n == 0:
        return np.array([])
    y_smooth = np.zeros(n, dtype=float)
    for i in range(n):
        d = np.abs(x - x[i])
        r = np.percentile(d, frac * 100.0)
        if r <= 0:
            y_smooth[i] = y[i]
            continue
        w = np.clip(1.0 - (d / r) ** 3, 0.0, 1.0) ** 3
        # Weighted linear regression
        b00 = np.sum(w)
        b01 = np.sum(w * x)
        b11 = np.sum(w * x * x)
        a0 = np.sum(w * y)
        a1 = np.sum(w * x * y)
        det = b00 * b11 - b01 * b01
        if abs(det) < 1e-12:
            y_smooth[i] = y[i]
        else:
            beta0 = (b11 * a0 - b01 * a1) / det
            beta1 = (b00 * a1 - b01 * a0) / det
            y_smooth[i] = beta0 + beta1 * x[i]
    return y_smooth


def _fit_spline_smooth(x: np.ndarray, y: np.ndarray, df: int = 10):
    """Fit a LOWESS smoother with adaptive bandwidth.

    Returns a callable that predicts at new x values using linear interpolation
    over the fitted curve (piecewise linear), which avoids the stepwise behavior
    of nearest-neighbor lookup.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return lambda x_new: np.zeros_like(np.asarray(x_new, dtype=float))

    x = x[mask]
    y = y[mask]

    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    # Adaptive fraction based on desired degrees of freedom
    frac = min(1.0, max(0.1, float(df) / float(len(x_sorted)))) if len(x_sorted) > 0 else 0.3
    y_fit_sorted = _fit_lowess_1d(x_sorted, y_sorted, frac=frac)

    # Collapse duplicate x values so `np.interp` is well-defined.
    x_unique, inv = np.unique(x_sorted, return_inverse=True)
    if len(x_unique) == 1:
        y_unique = np.array([float(np.nanmean(y_fit_sorted))], dtype=float)
    else:
        counts = np.bincount(inv)
        sums = np.bincount(inv, weights=y_fit_sorted)
        y_unique = sums / np.maximum(counts, 1)

    def predict(x_new: np.ndarray) -> np.ndarray:
        x_new = np.asarray(x_new, dtype=float)
        if x_new.ndim == 0:
            x_new = np.array([x_new])
        if x_unique.size == 0:
            return np.zeros_like(x_new, dtype=float)
        if x_unique.size == 1:
            return np.full_like(x_new, float(y_unique[0]), dtype=float)
        return np.interp(x_new, x_unique, y_unique)

    return predict


def _fit_linear_model(x: np.ndarray, y: np.ndarray):
    """Fit y ~ a + b*x via least squares and return predictor."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 2:
        return lambda x_new: np.zeros_like(np.asarray(x_new, dtype=float))
    x_m = x[mask]
    y_m = y[mask]
    x_mean = float(np.mean(x_m))
    y_mean = float(np.mean(y_m))
    denom = float(np.sum((x_m - x_mean) ** 2))
    if denom <= 1e-12:
        slope = 0.0
    else:
        slope = float(np.sum((x_m - x_mean) * (y_m - y_mean)) / denom)
    intercept = y_mean - slope * x_mean

    def predict(x_new: np.ndarray) -> np.ndarray:
        x_new = np.asarray(x_new, dtype=float)
        return intercept + slope * x_new

    return predict


def _fit_robust_linear_model(x: np.ndarray, y: np.ndarray, *, epsilon: float = 1.35):
    """Fit y ~ a + b*x using Huber regression and return predictor."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 2:
        return lambda x_new: np.zeros_like(np.asarray(x_new, dtype=float))

    x_m = x[mask].reshape(-1, 1)
    y_m = y[mask]
    try:
        from sklearn.linear_model import HuberRegressor

        model = HuberRegressor(epsilon=float(epsilon), alpha=0.0, fit_intercept=True, max_iter=200)
        model.fit(x_m, y_m)
        intercept = float(model.intercept_)
        slope = float(model.coef_[0]) if getattr(model, "coef_", None) is not None else 0.0
    except Exception:
        return _fit_linear_model(x, y)

    def predict(x_new: np.ndarray) -> np.ndarray:
        x_new = np.asarray(x_new, dtype=float)
        return intercept + slope * x_new

    return predict


def _fit_drift_models_core(
    tight_matches: pd.DataFrame,
    cfg: MassSightConfig,
) -> Dict[str, object]:
    """Fit drift correction models on tight matches."""
    if tight_matches.empty or len(tight_matches) < 10:
        # Return identity functions if insufficient data
        return {
            "rt_model": lambda x: np.zeros_like(np.asarray(x, dtype=float)),
            "rt_model_x": "rt",
            "ppm_model": lambda x: np.zeros_like(np.asarray(x, dtype=float)),
            "int_model": lambda x: np.zeros_like(np.asarray(x, dtype=float)),
            "roi_model": lambda x: np.zeros_like(np.asarray(x, dtype=float)),
        }

    rt_x_mode = str(getattr(cfg, "rt_drift_x", "rt") or "rt").lower().strip()
    if rt_x_mode in {"rt", "time"}:
        rt_x = tight_matches["rt_x"].to_numpy(dtype=float)
    elif rt_x_mode in {"ro", "roi", "rank"}:
        rt_x = tight_matches["roi_x"].to_numpy(dtype=float)
    else:
        raise ValueError(f"Unsupported rt_drift_x: {cfg.rt_drift_x!r}")

    delta_rt = tight_matches["delta_rt"].to_numpy(dtype=float)
    rt_mode = str(getattr(cfg, "rt_drift_model", "none") or "none").lower().strip()
    if rt_mode in {"none", "off", "identity"}:
        rt_model = lambda x: np.zeros_like(np.asarray(x, dtype=float))
    elif rt_mode in {"lowess", "smooth"}:
        rt_model = _fit_spline_smooth(rt_x, delta_rt, df=int(cfg.rt_spline_df))
    elif rt_mode in {"linear", "ols"}:
        rt_model = _fit_linear_model(rt_x, delta_rt)
    elif rt_mode in {"robust-linear", "robust_linear", "huber"}:
        rt_model = _fit_robust_linear_model(rt_x, delta_rt)
    else:
        raise ValueError(f"Unsupported rt_drift_model: {cfg.rt_drift_model!r}")

    mz_x = tight_matches["mz_x"].to_numpy(dtype=float)
    delta_ppm = tight_matches["delta_ppm"].to_numpy(dtype=float)
    ppm_mode = str(getattr(cfg, "ppm_drift_model", "lowess") or "lowess").lower().strip()
    if ppm_mode in {"none", "off", "identity"}:
        ppm_model = lambda x: np.zeros_like(np.asarray(x, dtype=float))
    elif ppm_mode in {"lowess", "smooth"}:
        ppm_model = _fit_spline_smooth(mz_x, delta_ppm, df=int(cfg.ppm_spline_df))
    elif ppm_mode in {"linear", "ols"}:
        ppm_model = _fit_linear_model(mz_x, delta_ppm)
    elif ppm_mode in {"robust-linear", "robust_linear", "huber"}:
        ppm_model = _fit_robust_linear_model(mz_x, delta_ppm)
    else:
        raise ValueError(f"Unsupported ppm_drift_model: {cfg.ppm_drift_model!r}")

    int_x = tight_matches["int_x"].to_numpy(dtype=float)
    int_diff = tight_matches["intensity_log_diff"].to_numpy(dtype=float)
    int_model = _fit_linear_model(int_x, int_diff)

    roi_x = tight_matches["roi_x"].to_numpy(dtype=float)
    delta_roi = tight_matches["delta_roi"].to_numpy(dtype=float)
    roi_mode = str(getattr(cfg, "roi_drift_model", "lowess") or "lowess").lower().strip()
    if roi_mode in {"none", "off", "identity"}:
        roi_model = lambda x: np.zeros_like(np.asarray(x, dtype=float))
    elif roi_mode in {"lowess", "smooth"}:
        roi_model = _fit_spline_smooth(roi_x, delta_roi, df=int(cfg.roi_spline_df))
    elif roi_mode in {"linear", "ols"}:
        roi_model = _fit_linear_model(roi_x, delta_roi)
    elif roi_mode in {"robust-linear", "robust_linear", "huber"}:
        roi_model = _fit_robust_linear_model(roi_x, delta_roi)
    else:
        raise ValueError(f"Unsupported roi_drift_model: {cfg.roi_drift_model!r}")

    return {
        "rt_model": rt_model,
        "rt_model_x": rt_x_mode,
        "ppm_model": ppm_model,
        "int_model": int_model,
        "roi_model": roi_model,
    }


def _bootstrap_resample_df(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    if df.empty:
        return df
    idx = rng.integers(0, len(df), size=len(df), dtype=np.int64)
    return df.iloc[idx].reset_index(drop=True)


def _fit_drift_models(
    tight_matches: pd.DataFrame,
    cfg: MassSightConfig,
) -> Dict[str, object]:
    """Fit drift correction models (with optional anchor bootstrap)."""
    bootstrap_fit = bool(getattr(cfg, "drift_fit_bootstrap", False))
    seed = int(getattr(cfg, "drift_bootstrap_seed", 0) or 0)

    if tight_matches.empty or len(tight_matches) < 10:
        models = _fit_drift_models_core(tight_matches, cfg)
        models.update(
            {
                "drift_fit_bootstrap": False,
                "drift_bootstrap_seed": seed,
            }
        )
        return models

    rng = np.random.default_rng(seed)
    tight_main = _bootstrap_resample_df(tight_matches, rng) if bootstrap_fit else tight_matches
    models = _fit_drift_models_core(tight_main, cfg)
    models.update(
        {
            "drift_fit_bootstrap": bool(bootstrap_fit),
            "drift_bootstrap_seed": int(seed),
        }
    )

    return models


def _apply_drift_correction(
    cand: pd.DataFrame,
    drift_models: Dict[str, object],
) -> pd.DataFrame:
    """Apply drift correction and compute residuals."""
    cand = cand.copy()

    rt_pred_x = cand["rt1"].to_numpy(dtype=float)
    rt_model_x = str(drift_models.get("rt_model_x", "rt") or "rt").lower().strip()
    if rt_model_x in {"ro", "roi", "rank"}:
        rt_pred_x = cand["roi_x"].to_numpy(dtype=float)
    rt_fit = drift_models["rt_model"](rt_pred_x)
    rt1 = cand["rt1"].to_numpy(dtype=float)
    delta_rt = cand["rt2"].to_numpy(dtype=float) - rt1
    cand["rt_fit"] = rt_fit
    cand["rt_resid"] = delta_rt - rt_fit

    mz_x = cand["mz1"].to_numpy(dtype=float)
    mz_y = cand["mz2"].to_numpy(dtype=float)
    if "mz_shift_da" in cand.columns:
        shift = cand["mz_shift_da"].to_numpy(dtype=float)
        shift = np.where(np.isfinite(shift), shift, 0.0)
    else:
        shift = np.zeros(len(cand), dtype=float)
    mz_x_eff = mz_x + shift
    cand["mz1_eff"] = mz_x_eff
    denom = np.maximum(0.5 * (mz_x_eff + mz_y), 1e-12)
    delta_ppm = (mz_y - mz_x_eff) / denom * 1e6
    ppm_fit = drift_models["ppm_model"](mz_x_eff)
    cand["delta_ppm_raw"] = delta_ppm
    cand["ppm_fit"] = ppm_fit
    cand["ppm_resid"] = delta_ppm - ppm_fit

    int_x = cand["int1"].to_numpy(dtype=float)
    int_fit = drift_models["int_model"](int_x)
    cand["int_fit"] = int_fit
    cand["int_resid"] = cand["intensity_log_diff"].to_numpy(dtype=float) - int_fit

    roi_x = cand["roi_x"].to_numpy(dtype=float)
    roi_fit = drift_models["roi_model"](roi_x)
    delta_roi = cand["delta_roi"].to_numpy(dtype=float)
    cand["roi_fit"] = roi_fit
    cand["roi_resid"] = delta_roi - roi_fit

    return cand


def _mad_scale(x: np.ndarray, floor: float) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float(max(floor, 1e-6))
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.std(x)) if np.std(x) > 0 else float(floor)
    return float(max(scale, floor))


def _student_t_logpdf(x: np.ndarray, df: float, scale: float) -> np.ndarray:
    if scale <= 0 or not np.isfinite(scale):
        scale = 1.0
    v = float(df)
    z = x / scale
    log_norm = math.lgamma((v + 1.0) / 2.0) - math.lgamma(v / 2.0) - 0.5 * math.log(v * math.pi) - math.log(scale)
    return log_norm - 0.5 * (v + 1.0) * np.log1p((z ** 2) / v)


def _anchor_residuals(
    tight: pd.DataFrame,
    drift_models: Dict[str, object],
) -> Dict[str, np.ndarray]:
    if tight.empty:
        return {}
    rt_model_x = str(drift_models.get("rt_model_x", "rt") or "rt").lower().strip()
    if rt_model_x in {"ro", "roi", "rank"}:
        rt_x = tight["roi_x"].to_numpy(dtype=float)
    else:
        rt_x = tight["rt_x"].to_numpy(dtype=float)
    mz_x = tight["mz_x"].to_numpy(dtype=float)
    roi_x = tight["roi_x"].to_numpy(dtype=float)
    int_x = tight["int_x"].to_numpy(dtype=float)

    rt_fit = drift_models["rt_model"](rt_x)
    ppm_fit = drift_models["ppm_model"](mz_x)
    int_fit = drift_models["int_model"](int_x)
    roi_fit = drift_models["roi_model"](roi_x)

    return {
        "ppm": tight["delta_ppm"].to_numpy(dtype=float) - ppm_fit,
        "rt": tight["delta_rt"].to_numpy(dtype=float) - rt_fit,
        "roi": tight["delta_roi"].to_numpy(dtype=float) - roi_fit,
        "int": tight["intensity_log_diff"].to_numpy(dtype=float) - int_fit,
    }


def _compute_anchor_scales(
    tight: pd.DataFrame,
    drift_models: Dict[str, object],
    cfg: MassSightConfig,
) -> Dict[str, float]:
    if tight.empty or len(tight) < int(cfg.anchor_scale_min_count):
        return {}
    resid = _anchor_residuals(tight, drift_models)
    if not resid:
        return {}
    scales = {
        "ppm": _mad_scale(resid["ppm"], cfg.scale_floor),
        "rt": _mad_scale(resid["rt"], cfg.scale_floor),
        "roi": _mad_scale(resid["roi"], cfg.scale_floor),
        "int": _mad_scale(resid["int"], cfg.scale_floor),
    }
    return scales


def _compute_loglik(
    cand: pd.DataFrame,
    cfg: MassSightConfig,
    *,
    scales_override: Optional[Dict[str, float]] = None,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    cand = cand.copy()
    scales: Dict[str, float] = {} if scales_override is None else dict(scales_override)

    ppm_resid = cand["ppm_resid"].to_numpy(dtype=float) if "ppm_resid" in cand.columns else np.zeros(len(cand))
    rt_resid = cand["rt_resid"].to_numpy(dtype=float) if "rt_resid" in cand.columns else np.zeros(len(cand))
    roi_resid = cand["roi_resid"].to_numpy(dtype=float) if "roi_resid" in cand.columns else np.zeros(len(cand))
    int_resid = cand["int_resid"].to_numpy(dtype=float) if "int_resid" in cand.columns else np.zeros(len(cand))

    if "ppm" not in scales:
        scales["ppm"] = _mad_scale(ppm_resid, cfg.scale_floor)
    if "rt" not in scales:
        scales["rt"] = _mad_scale(rt_resid, cfg.scale_floor)
    if "roi" not in scales:
        scales["roi"] = _mad_scale(roi_resid, cfg.scale_floor)
    if "int" not in scales:
        scales["int"] = _mad_scale(int_resid, cfg.scale_floor)

    ll = np.zeros(len(cand), dtype=float)

    if cfg.w_ppm != 0:
        ll += float(cfg.w_ppm) * _student_t_logpdf(ppm_resid, cfg.t_df, float(scales["ppm"]))
    if cfg.w_rt != 0:
        ll += float(cfg.w_rt) * _student_t_logpdf(rt_resid, cfg.t_df, float(scales["rt"]))
    if cfg.w_roi != 0:
        ll += float(cfg.w_roi) * _student_t_logpdf(roi_resid, cfg.t_df, float(scales["roi"]))
    if cfg.use_intensity and cfg.w_int != 0:
        ll += float(cfg.w_int) * _student_t_logpdf(int_resid, cfg.t_df, float(scales["int"]))

    shift_penalty = float(getattr(cfg, "mz_shift_penalty", 0.0) or 0.0)
    if shift_penalty != 0.0 and "mz_shift_da" in cand.columns:
        shift = cand["mz_shift_da"].to_numpy(dtype=float)
        ll = ll - shift_penalty * (np.abs(shift) > 1e-12).astype(float)

    cand["loglik"] = ll
    return cand, scales


def _sinkhorn_sparse(
    id1: np.ndarray,
    id2: np.ndarray,
    loglik: np.ndarray,
    n1: int,
    n2: int,
    cfg: MassSightConfig,
    *,
    a: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
    tau_a: float = 1.0,
    tau_b: float = 1.0,
) -> Tuple[np.ndarray, Dict[str, float]]:
    eps = float(cfg.ot_epsilon)
    min_mass = float(cfg.ot_min_mass)
    max_iter = int(cfg.ot_max_iter)
    tol = float(cfg.ot_tol)

    if eps <= 0 or not np.isfinite(eps):
        eps = 1.0

    loglik = loglik.copy()
    loglik[~np.isfinite(loglik)] = -np.inf

    # Row-wise stabilization
    row_max = pd.Series(loglik).groupby(id1).transform("max").to_numpy(dtype=float)
    row_max[~np.isfinite(row_max)] = 0.0
    logk = (loglik - row_max) / eps
    k = np.exp(logk)

    if a is None:
        row_has = np.bincount(id1, minlength=n1) > 0
        a = np.zeros(n1, dtype=float)
        if row_has.any():
            a[row_has] = 1.0 / float(row_has.sum())
    else:
        row_has = a > 0

    if b is None:
        col_has = np.bincount(id2, minlength=n2) > 0
        b = np.zeros(n2, dtype=float)
        if col_has.any():
            b[col_has] = 1.0 / float(col_has.sum())
    else:
        col_has = b > 0

    u = np.ones(n1, dtype=float)
    v = np.ones(n2, dtype=float)

    row_err = float("nan")
    col_err = float("nan")
    unbalanced = (abs(float(tau_a) - 1.0) > 1e-12) or (abs(float(tau_b) - 1.0) > 1e-12)
    tau_a = float(tau_a)
    tau_b = float(tau_b)
    if not np.isfinite(tau_a) or tau_a < 0:
        tau_a = 1.0
    if not np.isfinite(tau_b) or tau_b < 0:
        tau_b = 1.0

    for it in range(max_iter):
        u_prev = u.copy() if unbalanced else None
        v_prev = v.copy() if unbalanced else None

        r = np.bincount(id1, weights=k * v[id2], minlength=n1)
        u = np.zeros(n1, dtype=float)
        if unbalanced:
            u[row_has] = (a[row_has] / np.maximum(r[row_has], min_mass)) ** tau_a
        else:
            u[row_has] = a[row_has] / np.maximum(r[row_has], min_mass)

        c = np.bincount(id2, weights=k * u[id1], minlength=n2)
        v = np.zeros(n2, dtype=float)
        if unbalanced:
            v[col_has] = (b[col_has] / np.maximum(c[col_has], min_mass)) ** tau_b
        else:
            v[col_has] = b[col_has] / np.maximum(c[col_has], min_mass)

        if it % 10 == 0 or it == max_iter - 1:
            if unbalanced:
                du = float(np.nanmax(np.abs(u - (u_prev if u_prev is not None else u))))
                dv = float(np.nanmax(np.abs(v - (v_prev if v_prev is not None else v))))
                row_err = du
                col_err = dv
                if du <= tol and dv <= tol:
                    break
            else:
                r_check = np.bincount(id1, weights=k * v[id2], minlength=n1)
                row_sum = u * r_check
                c_check = np.bincount(id2, weights=k * u[id1], minlength=n2)
                col_sum = v * c_check
                if row_has.any():
                    row_err = float(np.nanmax(np.abs(row_sum[row_has] - a[row_has])))
                if col_has.any():
                    col_err = float(np.nanmax(np.abs(col_sum[col_has] - b[col_has])))
                if row_err <= tol and col_err <= tol:
                    break

    weights = k * u[id1] * v[id2]
    summary = {
        "ot_iter": float(it + 1),
        "ot_row_err": float(row_err),
        "ot_col_err": float(col_err),
        "ot_epsilon": float(eps),
        "ot_tau_row": float(tau_a),
        "ot_tau_col": float(tau_b),
        "ot_unbalanced": float(1.0 if unbalanced else 0.0),
    }
    return weights, summary


def _assignment_top1(cand: pd.DataFrame, score_col: str) -> pd.DataFrame:
    d = cand.sort_values(["id1", score_col], ascending=[True, False]).copy()
    top = d.groupby("id1").head(1)[["id1", "id2", score_col]].rename(columns={score_col: "prob"})
    return top.reset_index(drop=True)


def _solve_ot(
    cand: pd.DataFrame,
    loglik: np.ndarray,
    n1: int,
    n2: int,
    cfg: MassSightConfig,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Dict[str, float], Optional[float]]:
    id1 = cand["id1"].to_numpy(dtype=int)
    id2 = cand["id2"].to_numpy(dtype=int)
    loglik = np.asarray(loglik, dtype=float)

    null_mass = float(cfg.null_mass)
    null_mass = min(max(null_mass, 0.0), 0.9)
    a = None
    b = None
    null_loglik = None
    null_col = None
    n2_ext = n2
    n_edges_real = len(cand)
    row_has = np.bincount(id1, minlength=n1) > 0

    if cfg.allow_unmatched:
        col_has = np.bincount(id2, minlength=n2) > 0
        n_rows = int(row_has.sum())
        n_cols = int(col_has.sum())

        a = np.zeros(n1, dtype=float)
        if n_rows > 0:
            a[row_has] = 1.0 / float(n_rows)

        b_real = np.zeros(n2, dtype=float)
        if n_cols > 0:
            b_real[col_has] = (1.0 - null_mass) / float(n_cols)

        null_col = n2
        n2_ext = n2 + 1
        b = np.concatenate([b_real, np.array([null_mass], dtype=float)])

        if cfg.null_loglik is not None:
            null_loglik = float(cfg.null_loglik)
        else:
            finite = loglik[np.isfinite(loglik)]
            if finite.size:
                null_loglik = float(np.nanquantile(finite, float(cfg.null_loglik_quantile)))
            else:
                null_loglik = float(-10.0)

        null_rows = np.where(row_has)[0]
        id1 = np.concatenate([id1, null_rows.astype(int)], axis=0)
        id2 = np.concatenate([id2, np.full(len(null_rows), null_col, dtype=int)], axis=0)
        loglik = np.concatenate([loglik, np.full(len(null_rows), null_loglik, dtype=float)], axis=0)

    mode = str(getattr(cfg, "ot_mode", "balanced") or "balanced").strip().lower()
    mode = mode.replace("-", "_")
    if mode not in {"balanced", "unbalanced", "semi_relaxed"}:
        raise ValueError(f"Unsupported ot_mode: {cfg.ot_mode!r}")

    # --- Semi-relaxed OT: constrain rows only (columns have infinite capacity). ---
    # This yields a per-row softmax distribution at temperature `ot_epsilon` (plus null option).
    if mode == "semi_relaxed":
        # Row marginals used only to produce a transport mass (for diagnostics / compatibility).
        if a is None:
            row_has2 = np.bincount(id1[:n_edges_real], minlength=n1) > 0
            a = np.zeros(n1, dtype=float)
            if row_has2.any():
                a[row_has2] = 1.0 / float(row_has2.sum())

        # Probabilities over real edges (and per-row null via p0).
        p_row, p0 = _row_softmax_sparse(
            id1[:n_edges_real],
            loglik[:n_edges_real],
            n1,
            null_loglik=float(null_loglik) if (cfg.allow_unmatched and null_loglik is not None) else None,
            temperature=float(cfg.ot_epsilon),
        )
        weights_real = p_row * a[id1[:n_edges_real]]
        ot_summary = {
            "ot_mode": "semi_relaxed",
            "ot_iter": 1.0,
            "ot_row_err": float("nan"),
            "ot_col_err": float("nan"),
            "ot_epsilon": float(cfg.ot_epsilon),
            "ot_tau_row": float("nan"),
            "ot_tau_col": float("nan"),
            "ot_unbalanced": 0.0,
        }
        ot_summary.update({
            "allow_unmatched": float(1.0 if cfg.allow_unmatched else 0.0),
            "null_mass": float(null_mass if cfg.allow_unmatched else np.nan),
            "null_loglik": float(null_loglik) if cfg.allow_unmatched else float("nan"),
        })
        return weights_real, p_row, p0, ot_summary, null_loglik

    # --- Balanced / unbalanced OT via Sinkhorn on a sparse kernel. ---
    tau_a = 1.0
    tau_b = 1.0
    if mode == "unbalanced":
        lam_row = float(getattr(cfg, "ot_unbalanced_lambda_row", 10.0) or 10.0)
        lam_col = float(getattr(cfg, "ot_unbalanced_lambda_col", 1.0) or 1.0)
        eps = float(cfg.ot_epsilon)
        if not np.isfinite(lam_row) or lam_row < 0:
            lam_row = 10.0
        if not np.isfinite(lam_col) or lam_col < 0:
            lam_col = 1.0
        if not np.isfinite(eps) or eps <= 0:
            eps = 1.0
        tau_a = float(lam_row / (lam_row + eps))
        tau_b = float(lam_col / (lam_col + eps))

    weights, ot_summary = _sinkhorn_sparse(id1, id2, loglik, n1, n2_ext, cfg, a=a, b=b, tau_a=tau_a, tau_b=tau_b)

    weights_real = weights[:n_edges_real]
    row_real = np.bincount(id1[:n_edges_real], weights=weights_real, minlength=n1)
    row_null = np.zeros(n1, dtype=float)
    if cfg.allow_unmatched and null_col is not None:
        null_weights = weights[n_edges_real:]
        null_rows = id1[n_edges_real:]
        row_null = np.bincount(null_rows, weights=null_weights, minlength=n1)

    row_total = row_real + row_null
    denom = row_total[id1[:n_edges_real]]
    denom = np.where(denom > 0, denom, 1.0)
    p_row_ot = weights_real / denom
    p0 = None
    if cfg.allow_unmatched:
        p0 = np.zeros(n1, dtype=float)
        nonzero = row_total > 0
        p0[nonzero] = row_null[nonzero] / row_total[nonzero]

    ot_summary.update({
        "ot_mode": str(mode),
        "allow_unmatched": float(1.0 if cfg.allow_unmatched else 0.0),
        "null_mass": float(null_mass if cfg.allow_unmatched else np.nan),
        "null_loglik": float(null_loglik) if cfg.allow_unmatched else float("nan"),
    })

    return weights_real, p_row_ot, p0, ot_summary, null_loglik


def _row_softmax_sparse(
    id1: np.ndarray,
    score: np.ndarray,
    n1: int,
    *,
    null_loglik: Optional[float] = None,
    temperature: float = 1.0,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Row-wise softmax over a sparse edge list, optionally including a null option.

    Returns:
      - p_row: per-edge probabilities, normalized within each id1 (and including null if provided)
      - p0: per-row null probability (length n1) if null_loglik is not None
    """
    id1 = np.asarray(id1, dtype=int)
    score = np.asarray(score, dtype=float)
    n1 = int(n1)
    if n1 <= 0 or score.size == 0:
        p0 = None
        if null_loglik is not None and n1 > 0:
            p0 = np.ones(n1, dtype=float)
        return np.zeros_like(score, dtype=float), p0

    temp = float(temperature)
    if not np.isfinite(temp) or temp <= 0:
        temp = 1.0

    s = score.copy()
    s[~np.isfinite(s)] = -np.inf

    # Stable row max via numpy (avoid pandas overhead).
    row_max = np.full(n1, -np.inf, dtype=float)
    np.maximum.at(row_max, id1, s)
    row_max[~np.isfinite(row_max)] = 0.0

    exp_scores = np.exp((s - row_max[id1]) / temp)
    sum_exp = np.bincount(id1, weights=exp_scores, minlength=n1).astype(float)

    if null_loglik is None:
        denom = sum_exp[id1]
        denom = np.where(denom > 0, denom, 1.0)
        return exp_scores / denom, None

    nl = float(null_loglik)
    # Include the null option for every row (even those with no edges).
    exp_null = np.exp((nl - row_max) / temp)
    total = sum_exp + exp_null
    total_safe = np.where(total > 0, total, 1.0)
    p0 = exp_null / total_safe
    return exp_scores / total_safe[id1], p0


def _row_entropy_sparse(id1: np.ndarray, p_row: np.ndarray, n1: int, *, p0: Optional[np.ndarray] = None) -> np.ndarray:
    """Per-row Shannon entropy for a sparse row distribution (optionally with null mass)."""
    id1 = np.asarray(id1, dtype=int)
    p_row = np.asarray(p_row, dtype=float)
    n1 = int(n1)
    if n1 <= 0:
        return np.zeros(0, dtype=float)
    ent = np.zeros(n1, dtype=float)
    ok = np.isfinite(p_row) & (p_row > 0)
    if np.any(ok):
        term = p_row[ok] * np.log(p_row[ok])
        np.add.at(ent, id1[ok], term)
    if p0 is not None:
        p0 = np.asarray(p0, dtype=float)
        ok0 = np.isfinite(p0) & (p0 > 0)
        ent[ok0] += p0[ok0] * np.log(p0[ok0])
    return -ent


def _build_row_maps(cand: pd.DataFrame, n1: int, prob_col: str) -> List[Dict[int, float]]:
    maps: List[Dict[int, float]] = [dict() for _ in range(n1)]
    for row in cand.itertuples(index=False):
        prob = getattr(row, prob_col)
        if prob <= 0:
            continue
        maps[int(row.id1)][int(row.id2)] = float(prob)
    return maps


def _compute_structure_bonus(
    cand: pd.DataFrame,
    row_maps: List[Dict[int, float]],
    graph1: StructureGraph,
    graph2: StructureGraph,
    cfg: MassSightConfig,
) -> np.ndarray:
    mode = str(cfg.structure_bonus_mode).lower().strip()
    if mode not in {"mass", "logmass", "gw"}:
        raise ValueError(f"Unsupported structure_bonus_mode: {cfg.structure_bonus_mode!r}")

    id1 = cand["id1"].to_numpy(dtype=int)
    id2 = cand["id2"].to_numpy(dtype=int)
    bonus = np.zeros(len(cand), dtype=float)

    # Fast path: rewrite the inner loop to iterate over each neighbor's sparse
    # candidate map rather than over all (k^2) neighbor pairs.
    n1 = graph1.neighbors
    w1 = graph1.weights
    n2 = graph2.neighbors
    w2 = graph2.weights
    if mode == "gw":
        sumsq1 = np.sum(np.asarray(w1, dtype=float) ** 2, axis=1)
        sumsq2 = np.sum(np.asarray(w2, dtype=float) ** 2, axis=1)
        # Normalize the scale of the GW term so that structure_alpha is roughly
        # comparable across different k / graph densities.
        # Typical sumsq ~ 1/k for normalized kNN weights; invert to get O(1) scores.
        typical = 0.5 * (float(np.nanmedian(sumsq1)) + float(np.nanmedian(sumsq2)))
        gw_scale = 1.0 / max(typical, float(cfg.structure_eps))

    # Precompute neighbor-weight lookup for graph2 so we can score membership in O(1).
    # This enables using larger k without quadratic slowdowns.
    nbr_weight2: List[Dict[int, float]] = []
    for nbrs, wts in zip(n2, w2):
        m: Dict[int, float] = {}
        for jj, ww in zip(nbrs, wts):
            if jj < 0 or ww == 0.0:
                continue
            m[int(jj)] = float(ww)
        nbr_weight2.append(m)

    for idx, (i, j) in enumerate(zip(id1, id2)):
        cross = 0.0
        wmap2 = nbr_weight2[int(j)]
        for i2, wi in zip(n1[int(i)], w1[int(i)]):
            if i2 < 0 or wi == 0.0:
                continue
            row_map = row_maps[int(i2)]
            if not row_map:
                continue
            inner = 0.0
            for j2, p in row_map.items():
                wj = wmap2.get(int(j2))
                if wj:
                    inner += wj * float(p)
            cross += float(wi) * inner
        if mode == "gw":
            # Approximate squared-loss GW gradient (up to row-constant terms):
            # cost(i,j) = sum_k C1(i,k)^2 + sum_l C2(j,l)^2 - 2 * <C1(i,*), Π C2(j,*)>
            cost = float(sumsq1[int(i)] + sumsq2[int(j)] - 2.0 * cross)
            bonus[idx] = -gw_scale * cost
        else:
            bonus[idx] = cross

    if mode == "logmass":
        bonus = np.log(bonus + float(cfg.structure_eps))
    return bonus


def _prepare_match(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    Dict[str, object],
    Dict[str, float],
    Dict[str, float],
    Dict[str, object],
]:
    """Prepare massSight candidates with drift correction and local log-likelihood."""
    polarity = str(cfg.polarity).lower().strip()
    if polarity not in {"positive", "negative"}:
        raise ValueError(f"Unsupported polarity '{cfg.polarity}'.")

    cfg_eff = cfg
    rt_gate_mode = str(getattr(cfg, "rt_gate_mode", "none") or "none").lower().strip()
    if rt_gate_mode not in {"none", "adaptive-disable"}:
        raise ValueError(
            f"Unsupported rt_gate_mode: {cfg.rt_gate_mode!r} (expected 'none' or 'adaptive-disable')."
        )
    rt_gate_diag: Dict[str, object] = {
        "rt_gate_mode": rt_gate_mode,
        "rt_gate_rt_requested": float(cfg.rt) if cfg.rt is not None else 0.0,
        "rt_gate_disable_p95": float(cfg.rt_gate_disable_p95),
        "rt_gate_anchor_ppm": float(cfg.rt_gate_anchor_ppm)
        if cfg.rt_gate_anchor_ppm is not None and np.isfinite(float(cfg.rt_gate_anchor_ppm))
        else None,
        "rt_gate_min_count": int(cfg.rt_gate_min_count),
    }

    schema_cfg = SchemaConfig(
        mz_col=cfg.mz_col,
        rt_col=cfg.rt_col,
        intensity_col=cfg.intensity_col,
        annotation_col=cfg.annotation_col,
        compound_id_col=cfg.compound_id_col,
    )

    # Allow running in m/z-only mode even when RT is missing as a column (public software robustness).
    retention_mode_req = str(getattr(cfg_eff, "retention_mode", "auto") or "auto").lower().strip()
    if retention_mode_req in {"mz-only", "mz", "mz_only"}:
        retention_mode_req = "mz_only"
    if retention_mode_req in {"auto", "mz_only"}:
        rt_src_a = cfg.rt_col_study_a or schema_cfg.rt_col
        rt_src_b = cfg.rt_col_study_b or schema_cfg.rt_col
        if rt_src_a not in study_a.columns:
            study_a = study_a.copy()
            study_a[rt_src_a] = np.nan
        if rt_src_b not in study_b.columns:
            study_b = study_b.copy()
            study_b[rt_src_b] = np.nan

    study_a = normalize_schema(
        study_a,
        schema_cfg,
        mz_override=cfg.mz_col_study_a,
        rt_override=cfg.rt_col_study_a,
        annotation_override=cfg.annotation_col_study_a,
        compound_override=cfg.compound_id_col_study_a,
    )
    study_b = normalize_schema(
        study_b,
        schema_cfg,
        mz_override=cfg.mz_col_study_b,
        rt_override=cfg.rt_col_study_b,
        annotation_override=cfg.annotation_col_study_b,
        compound_override=cfg.compound_id_col_study_b,
    )

    # --- RT unit normalization (seconds -> minutes) ---
    if bool(getattr(cfg_eff, "normalize_rt_units", True)):
        scale_a, info_a = infer_rt_unit_scale_to_minutes(study_a["RT"].to_numpy(dtype=float))
        scale_b, info_b = infer_rt_unit_scale_to_minutes(study_b["RT"].to_numpy(dtype=float))
        rt_gate_diag["rt_unit_scale_study_a"] = float(scale_a)
        rt_gate_diag["rt_unit_scale_study_b"] = float(scale_b)
        rt_gate_diag["rt_unit_note_study_a"] = str(info_a.get("note", ""))
        rt_gate_diag["rt_unit_note_study_b"] = str(info_b.get("note", ""))
        if float(scale_a) != 1.0:
            study_a["RT"] = study_a["RT"].to_numpy(dtype=float) * float(scale_a)
        if float(scale_b) != 1.0:
            study_b["RT"] = study_b["RT"].to_numpy(dtype=float) * float(scale_b)

    cfg_eff, retention_diag = _maybe_disable_retention(study_a, study_b, cfg_eff)
    rt_gate_diag.update(retention_diag)

    # --- Per-pair RT transferability (fit linear transform or ignore retention) ---
    rt_transfer_mode = str(getattr(cfg_eff, "rt_transferability", "auto") or "auto").lower().strip()
    if rt_transfer_mode not in {"off", "auto"}:
        raise ValueError(
            f"Unsupported rt_transferability: {getattr(cfg_eff, 'rt_transferability', None)!r} "
            "(expected 'off' or 'auto')."
        )
    rt_gate_diag["rt_transferability"] = rt_transfer_mode
    rt_gate_diag["rt_policy"] = "off" if rt_transfer_mode == "off" else "unknown"
    rt_gate_diag["rt_scale"] = 1.0
    rt_gate_diag["rt_offset"] = 0.0
    rt_gate_diag["rt_note"] = ""

    retention_mode_req = str(retention_diag.get("retention_mode_requested", "auto") or "auto").lower().strip()
    retention_mode_eff = str(retention_diag.get("retention_mode_effective", retention_mode_req) or retention_mode_req)

    if (
        rt_transfer_mode == "auto"
        and retention_mode_eff != "mz_only"
        and (float(cfg_eff.w_rt) != 0.0 or float(cfg_eff.w_roi) != 0.0)
    ):
        rt_pair = _infer_pair_rt_transferability(study_a, study_b, cfg_eff)
        # Flatten under stable keys for downstream logging.
        for k, v in rt_pair.items():
            rt_gate_diag[str(k)] = v

        if str(rt_pair.get("rt_policy")) == "linear":
            scale = float(rt_pair.get("rt_scale", 1.0) or 1.0)
            offset = float(rt_pair.get("rt_offset", 0.0) or 0.0)
            study_b["RT"] = study_b["RT"].to_numpy(dtype=float) * scale + offset
            rt_gate_diag["rt_applied"] = True
        elif str(rt_pair.get("rt_policy")) == "ignore" and retention_mode_req == "auto":
            # In auto mode, treat non-transferable retention as uninformative and disable RT/ROI terms.
            cfg_eff = replace(
                cfg_eff,
                w_rt=0.0,
                w_roi=0.0,
                rt=0.0,
                tight_rt=0.0,
                tight_roi=0.0,
                rt_drift_model="none",
                roi_drift_model="none",
            )
            rt_gate_diag["retention_mode_effective"] = "mz_only"
            rt_gate_diag["retention_disable_reason"] = "rt_nontransferable_pair"

    if rt_gate_mode == "adaptive-disable" and float(cfg_eff.rt or 0.0) > 0:
        anchor_ppm = cfg_eff.rt_gate_anchor_ppm
        if anchor_ppm is None or not np.isfinite(float(anchor_ppm)) or float(anchor_ppm) <= 0:
            anchor_ppm = float(cfg_eff.tight_ppm)
        p95, n_anchors = _estimate_rt_drift_p95_mz_tight(
            study_a,
            study_b,
            ppm=float(anchor_ppm),
            min_count=int(cfg_eff.rt_gate_min_count),
        )
        rt_gate_diag["rt_gate_anchor_ppm"] = float(anchor_ppm)
        rt_gate_diag["rt_gate_p95_est"] = float(p95) if np.isfinite(p95) else None
        rt_gate_diag["rt_gate_n_anchors"] = int(n_anchors)
        disable = bool(np.isfinite(p95) and float(p95) > float(cfg_eff.rt_gate_disable_p95))
        rt_gate_diag["rt_gate_disabled"] = bool(disable)
        if disable:
            cfg_eff = replace(cfg_eff, rt=0.0)
        rt_gate_diag["rt_gate_rt_final"] = float(cfg_eff.rt) if cfg_eff.rt is not None else 0.0
    else:
        rt_gate_diag["rt_gate_p95_est"] = None
        rt_gate_diag["rt_gate_n_anchors"] = None
        rt_gate_diag["rt_gate_disabled"] = False
        rt_gate_diag["rt_gate_rt_final"] = float(cfg_eff.rt) if cfg_eff.rt is not None else 0.0

    shifts_da, mz_shift_diag = _select_mz_shifts(study_a, study_b, cfg_eff)
    mz_shift_diag.update(rt_gate_diag)
    cand = _build_candidate_pairs_shifted(study_a, study_b, cfg_eff, shifts_da)
    tight = _build_tight_matches_shifted(study_a, study_b, cfg_eff, shifts_da)
    if cand.empty:
        return study_a, study_b, cand, pd.DataFrame(), {}, {}, {}, mz_shift_diag

    drift_models = _fit_drift_models(tight, cfg_eff)
    cand = _apply_drift_correction(cand, drift_models)

    anchor_scales = _compute_anchor_scales(tight, drift_models, cfg_eff) if cfg_eff.use_anchor_scales else {}
    cand, scales = _compute_loglik(cand, cfg_eff, scales_override=anchor_scales or None)

    return study_a, study_b, cand, tight, drift_models, scales, anchor_scales, mz_shift_diag


def match_features(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    config: Optional[MassSightConfig] = None,
    *,
    expr_a: Optional[object] = None,
    expr_b: Optional[object] = None,
) -> MatchResult:
    cfg = config or MassSightConfig()

    study_a, study_b, cand, tight, drift_models, scales, anchor_scales, mz_shift_diag = _prepare_match(
        study_a, study_b, cfg
    )
    if cand.empty:
        empty = pd.DataFrame(
            {
                "id1": np.arange(len(study_a), dtype=int),
                "id2": -1,
                "decision": "no_candidate",
                "support": 1.0,
                "prob_raw": 0.0,
                "prob_match_raw": 0.0,
                "p_null": 0.0,
                "margin": 0.0,
            }
        )
        return MatchResult(
            candidates=cand,
            top1=empty,
            ot_summary={"n_candidates": 0.0},
            drift_params=mz_shift_diag,
        )

    cand = cand.copy()
    cand["loglik_local"] = cand["loglik"]
    cand["struct_bonus"] = 0.0
    cand["loglik_total"] = cand["loglik_local"]
    cand["loglik"] = cand["loglik_total"]

    n1 = len(study_a)
    n2 = len(study_b)

    cfg_use = cfg
    # Optional auto-structure fallback: when RT is deemed non-transferable across the pair,
    # within-study correlation structure can provide the missing disambiguation signal.
    try:
        rt_policy = str(mz_shift_diag.get("rt_policy", "") or "")
    except Exception:
        rt_policy = ""
    auto_alpha = float(getattr(cfg, "structure_alpha_on_nontransferable_rt", 0.0) or 0.0)
    if (
        (not bool(cfg.use_structure))
        and auto_alpha != 0.0
        and rt_policy == "ignore"
        and expr_a is not None
        and expr_b is not None
    ):
        cfg_use = replace(cfg, use_structure=True, structure_alpha=float(auto_alpha), structure_graph="local_rt")

    use_structure = bool(cfg_use.use_structure) and float(cfg_use.structure_alpha) != 0.0 and int(cfg_use.structure_iters) > 0
    graph1 = graph2 = None
    if use_structure:
        if expr_a is None or expr_b is None:
            raise ValueError("expr_a/expr_b are required when use_structure=True.")
        x1 = coerce_expression(expr_a, n1, "expr_a")
        x2 = coerce_expression(expr_b, n2, "expr_b")
        if x1.shape[0] >= int(cfg_use.structure_min_samples) and x2.shape[0] >= int(cfg_use.structure_min_samples):
            graph_mode = str(cfg_use.structure_graph).lower().strip()
            if graph_mode not in {"knn", "stable", "local_rt"}:
                raise ValueError(
                    f"Unsupported structure_graph: {cfg_use.structure_graph!r} (expected 'knn', 'stable', or 'local_rt')."
                )

            if graph_mode == "knn":
                graph1 = build_structure_graph(
                    x1,
                    k=int(cfg_use.structure_k),
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    eps=float(cfg_use.structure_eps),
                )
                graph2 = build_structure_graph(
                    x2,
                    k=int(cfg_use.structure_k),
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    eps=float(cfg_use.structure_eps),
                )
            elif graph_mode == "stable":
                graph1 = build_structure_graph_stable(
                    x1,
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    k_max=int(cfg_use.structure_k),
                    n_boot=int(cfg_use.structure_n_boot),
                    subsample_frac=float(cfg_use.structure_subsample_frac),
                    min_freq=float(cfg_use.structure_min_freq),
                    eps=float(cfg_use.structure_eps),
                    seed=0,
                )
                graph2 = build_structure_graph_stable(
                    x2,
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    k_max=int(cfg_use.structure_k),
                    n_boot=int(cfg_use.structure_n_boot),
                    subsample_frac=float(cfg_use.structure_subsample_frac),
                    min_freq=float(cfg_use.structure_min_freq),
                    eps=float(cfg_use.structure_eps),
                    seed=1,
                )
            else:
                rt1 = study_a["RT"].to_numpy(dtype=float, copy=False)
                rt2 = study_b["RT"].to_numpy(dtype=float, copy=False)
                graph1 = build_structure_graph_local_rt(
                    x1,
                    rt1,
                    k=int(cfg_use.structure_k),
                    rt_window=float(cfg_use.structure_rt_window),
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    max_candidates=int(cfg_use.structure_max_candidates),
                    eps=float(cfg_use.structure_eps),
                )
                graph2 = build_structure_graph_local_rt(
                    x2,
                    rt2,
                    k=int(cfg_use.structure_k),
                    rt_window=float(cfg_use.structure_rt_window),
                    corr_method=str(cfg_use.structure_corr).lower(),
                    corr_estimator=str(cfg_use.structure_corr_estimator).lower(),
                    use_abs=bool(cfg_use.structure_abs),
                    max_candidates=int(cfg_use.structure_max_candidates),
                    eps=float(cfg_use.structure_eps),
                )
        else:
            use_structure = False

    n_iters = int(cfg_use.structure_iters) if use_structure else 1
    ot_summary = {}
    p0 = None
    null_loglik = None
    for iter_idx in range(max(n_iters, 1)):
        cfg_ot = cfg_use
        if use_structure and cfg_use.structure_anneal and n_iters > 1:
            t = float(iter_idx) / float(max(n_iters - 1, 1))
            mult_start = float(cfg_use.structure_ot_epsilon_mult_start)
            mult_start = max(mult_start, 1.0)
            eps_iter = float(cfg_use.ot_epsilon) * (mult_start * (1.0 - t) + 1.0 * t)
            cfg_ot = replace(cfg_use, ot_epsilon=eps_iter)

        weights_real, p_row_ot, p0, ot_summary, null_loglik = _solve_ot(
            cand,
            cand["loglik_total"].to_numpy(dtype=float),
            n1,
            n2,
            cfg_ot,
        )
        cand["ot_weight"] = weights_real
        cand["p_row_ot"] = p_row_ot
        if cfg_use.allow_unmatched and p0 is not None:
            cand["p0_ot"] = cand["id1"].map(pd.Series(p0))

        if not use_structure or iter_idx == n_iters - 1 or graph1 is None or graph2 is None:
            break

        row_maps = _build_row_maps(cand, n1, "p_row_ot")
        bonus = _compute_structure_bonus(cand, row_maps, graph1, graph2, cfg_use)
        cand["struct_bonus"] = bonus
        alpha = float(cfg_use.structure_alpha)
        if cfg_use.structure_anneal and cfg_use.structure_alpha_ramp and n_iters > 1:
            alpha = alpha * float(iter_idx + 1) / float(n_iters - 1)

        if cfg_use.structure_adaptive_alpha and float(cfg_use.structure_margin_tau) > 0.0:
            ranked = cand.sort_values(["id1", "p_row_ot"], ascending=[True, False]).groupby("id1").head(2).copy()
            ranked["rank"] = ranked.groupby("id1").cumcount()
            pivot = ranked.pivot(index="id1", columns="rank", values="p_row_ot")
            top1 = pivot.get(0)
            top2 = pivot.get(1)
            if top1 is None:
                alpha_row = pd.Series(0.0, index=pd.RangeIndex(n1))
            else:
                margin = top1.fillna(0.0) - (top2.fillna(0.0) if top2 is not None else 0.0)
                tau = float(cfg_use.structure_margin_tau)
                # Smoothly downweight structure on confident rows (large top1-top2 margin)
                # without fully disabling it (important when structure provides global consistency).
                margin_v = margin.to_numpy(dtype=float)
                scale = 1.0 / (1.0 + (np.maximum(margin_v, 0.0) / tau))
                gamma = float(cfg_use.structure_margin_gamma)
                if gamma != 1.0:
                    scale = scale ** gamma
                alpha_row = pd.Series(scale, index=margin.index)
            alpha_edge = cand["id1"].map(alpha_row).fillna(0.0).to_numpy(dtype=float) * alpha
            cand["loglik_total"] = cand["loglik_local"] + alpha_edge * bonus
        else:
            cand["loglik_total"] = cand["loglik_local"] + alpha * bonus
        cand["loglik"] = cand["loglik_total"]

    # === Output ===
    # We always solve OT (above) to obtain a globally consistent soft coupling and to compute
    # a principled null utility. We report OT top‑1 (soft) with an explicit no‑match option.
    #
    # In addition, we expose a "local" correspondence distribution that does not impose
    # global 1‑1 competition: a row-wise softmax over `loglik_total` (plus the same null option).
    id1 = cand["id1"].to_numpy(dtype=int)
    p_row_local, p0_local = _row_softmax_sparse(
        id1,
        cand["loglik_total"].to_numpy(dtype=float),
        n1,
        null_loglik=(float(null_loglik) if (cfg_use.allow_unmatched and null_loglik is not None) else None),
    )
    cand["p_row_local"] = p_row_local
    if cfg_use.allow_unmatched and p0_local is not None:
        cand["p0_local"] = cand["id1"].map(pd.Series(p0_local))

    # Per-row uncertainty diagnostics (entropy + effective candidate count).
    row_ent_local = _row_entropy_sparse(id1, p_row_local, n1, p0=p0_local if cfg_use.allow_unmatched else None)
    cand["row_entropy_local"] = cand["id1"].map(pd.Series(row_ent_local))
    cand["row_neff_local"] = cand["id1"].map(pd.Series(np.exp(row_ent_local)))
    if "p_row_ot" in cand.columns:
        p_row_ot = cand["p_row_ot"].to_numpy(dtype=float)
        row_ent_ot = _row_entropy_sparse(id1, p_row_ot, n1, p0=p0 if cfg_use.allow_unmatched else None)
        cand["row_entropy_ot"] = cand["id1"].map(pd.Series(row_ent_ot))
        cand["row_neff_ot"] = cand["id1"].map(pd.Series(np.exp(row_ent_ot)))

    # Build per-row top‑1 and OT ambiguity margin (p1-p2).
    ranked = (
        cand.sort_values(["id1", "p_row_ot"], ascending=[True, False])
        .groupby("id1")
        .head(2)
        .copy()
    )
    ranked["rank"] = ranked.groupby("id1").cumcount()
    pivot = ranked.pivot(index="id1", columns="rank", values="p_row_ot")
    p1 = pivot.get(0, pd.Series(dtype=float)).reindex(range(int(n1))).fillna(0.0).to_numpy(dtype=float)
    p2 = pivot.get(1, pd.Series(dtype=float)).reindex(range(int(n1))).fillna(0.0).to_numpy(dtype=float)
    margin = np.maximum(p1 - p2, 0.0)

    best = ranked[ranked["rank"] == 0].loc[:, ["id1", "id2", "p_row_ot"]].copy()
    best = best.rename(columns={"id2": "id2_raw", "p_row_ot": "prob_match_raw"})
    top1 = pd.DataFrame({"id1": np.arange(int(n1), dtype=int)})
    top1 = top1.merge(best, how="left", on="id1")
    top1["id2_raw"] = top1["id2_raw"].fillna(-1).astype(int)
    top1["prob_match_raw"] = top1["prob_match_raw"].fillna(0.0)
    top1["margin"] = margin.astype(float)

    top1["id2"] = top1["id2_raw"].to_numpy(dtype=int)
    top1["decision"] = np.where(top1["id2"].to_numpy(dtype=int) >= 0, "match", "no_candidate")

    # OT null decision.
    if cfg_use.allow_unmatched and p0 is not None:
        p0_map = pd.Series(p0)
        top1["p_null"] = top1["id1"].map(p0_map).fillna(0.0)
        use_null = (top1["id2_raw"].to_numpy(dtype=int) >= 0) & (top1["p_null"] >= top1["prob_match_raw"])
        top1.loc[use_null, "id2"] = -1
        top1.loc[use_null, "decision"] = "null_ot"
    else:
        top1["p_null"] = 0.0

    # Final output probability: for matches report OT p1; for no-match report OT p_null.
    top1["prob_raw"] = top1["prob_match_raw"].to_numpy(dtype=float)
    no_match = top1["id2"].to_numpy(dtype=int) < 0
    if "p_null" in top1.columns:
        top1.loc[no_match, "prob_raw"] = top1.loc[no_match, "p_null"]

    if cfg_use.max_p0 is not None:
        top1 = top1[top1["p_null"] < float(cfg_use.max_p0)].copy()

    ot_summary.update({
        "n_candidates": float(len(cand)),
        "scale_ppm": float(scales.get("ppm", np.nan)),
        "scale_rt": float(scales.get("rt", np.nan)),
        "scale_roi": float(scales.get("roi", np.nan)),
        "scale_int": float(scales.get("int", np.nan)),
        "t_df": float(cfg_use.t_df),
        "anchor_scale_used": float(1.0 if anchor_scales else 0.0),
        "structure_enabled": float(1.0 if use_structure else 0.0),
        "structure_alpha": float(cfg_use.structure_alpha),
        "structure_iters": float(n_iters),
    })

    drift_params = {
        "tight_matches": len(tight),
        "drift_models": drift_models,
    }
    drift_params.update(mz_shift_diag)

    return MatchResult(
        candidates=cand,
        top1=top1,
        ot_summary=ot_summary,
        drift_params=drift_params,
    )


__all__ = [
    "MassSightConfig",
    "MatchResult",
    "match_features",
]
