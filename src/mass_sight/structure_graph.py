from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Tuple

import numpy as np

try:
    from scipy.stats import rankdata
except Exception:  # pragma: no cover
    rankdata = None  # type: ignore[assignment]

CorrMethod = Literal["pearson", "spearman"]
CorrEstimator = Literal["sample", "ledoit_wolf"]


@dataclass(frozen=True)
class StructureGraph:
    neighbors: np.ndarray  # shape: (n_features, k), -1 for padding
    weights: np.ndarray  # shape: (n_features, k), normalized (sum abs or sum)


def coerce_expression(expr: object, n_features: int, name: str) -> np.ndarray:
    """Return a 2D float array with shape (n_samples, n_features)."""
    if expr is None:
        raise ValueError(f"{name} is required when use_structure=True.")

    if hasattr(expr, "to_numpy"):
        arr = expr.to_numpy()
    else:
        arr = np.asarray(expr)

    if arr.ndim != 2:
        raise ValueError(f"{name} must be 2D, got shape {arr.shape}.")

    if arr.shape[1] == n_features:
        return arr.astype(float, copy=False)

    if arr.shape[0] == n_features and arr.shape[1] != n_features:
        return arr.T.astype(float, copy=False)

    raise ValueError(
        f"{name} has shape {arr.shape}, expected (n_samples, {n_features}) "
        f"or ({n_features}, n_samples)."
    )


def _spearman_rank(x: np.ndarray) -> np.ndarray:
    if rankdata is None:
        raise ImportError("scipy is required for spearman correlations.")
    try:
        return rankdata(x, axis=0)
    except TypeError:
        return np.apply_along_axis(rankdata, 0, x)


def build_structure_graph(
    expr: np.ndarray,
    *,
    k: int,
    corr_method: CorrMethod = "spearman",
    corr_estimator: CorrEstimator | str = "sample",
    use_abs: bool = True,
    eps: float = 1e-8,
) -> StructureGraph:
    """Build a sparse kNN correlation graph from samples×features matrix."""
    if expr.ndim != 2:
        raise ValueError(f"expr must be 2D, got shape {expr.shape}.")

    n_samples, n_features = expr.shape
    if n_features <= 1 or k <= 0:
        return StructureGraph(
            neighbors=np.full((n_features, 0), -1, dtype=int),
            weights=np.zeros((n_features, 0), dtype=float),
        )

    x = expr.astype(np.float32, copy=False)
    if corr_method == "spearman":
        x = _spearman_rank(x).astype(np.float32, copy=False)

    x = x - x.mean(axis=0, keepdims=True)

    est = str(corr_estimator).strip().lower().replace("-", "_")
    if est in {"sample", "empirical", "none"}:
        std = x.std(axis=0, ddof=1)
        std = np.where(std > eps, std, 1.0)
        x = x / std
        denom = float(max(n_samples - 1, 1))
        corr = (x.T @ x) / denom
    elif est in {"ledoit_wolf", "ledoitwolf", "lw"}:
        try:
            from sklearn.covariance import LedoitWolf  # type: ignore
        except Exception as e:  # pragma: no cover
            raise ImportError("scikit-learn is required for corr_estimator='ledoit_wolf'.") from e

        # Estimate a well-conditioned covariance, then convert to correlation.
        # This is especially useful when p >> n, where sample correlations are
        # unstable and produce noisy neighbor graphs.
        cov = LedoitWolf(assume_centered=True).fit(x).covariance_
        cov = np.asarray(cov, dtype=np.float32)
        d = np.sqrt(np.clip(np.diag(cov), eps, None))
        corr = cov / np.outer(d, d)
    else:
        raise ValueError(f"Unsupported corr_estimator: {corr_estimator!r} (expected 'sample' or 'ledoit_wolf').")

    score = np.abs(corr)
    np.fill_diagonal(score, -np.inf)

    k_eff = min(int(k), n_features - 1)
    idx = np.argpartition(score, -k_eff, axis=1)[:, -k_eff:]
    row_score = np.take_along_axis(score, idx, axis=1)
    order = np.argsort(row_score, axis=1)[:, ::-1]
    idx = np.take_along_axis(idx, order, axis=1)

    weights = np.take_along_axis(corr, idx, axis=1)
    if use_abs:
        weights = np.abs(weights)
        denom = weights.sum(axis=1, keepdims=True)
    else:
        denom = np.sum(np.abs(weights), axis=1, keepdims=True)
    denom = np.where(denom > eps, denom, 1.0)
    weights = weights / denom

    neighbors = idx.astype(int, copy=False)
    weights = weights.astype(float, copy=False)

    return StructureGraph(neighbors=neighbors, weights=weights)


def build_structure_graph_stable(
    expr: np.ndarray,
    *,
    corr_method: CorrMethod = "spearman",
    corr_estimator: CorrEstimator | str = "sample",
    use_abs: bool = True,
    k_max: int = 50,
    n_boot: int = 10,
    subsample_frac: float = 0.8,
    min_freq: float = 0.6,
    eps: float = 1e-8,
    seed: int = 0,
) -> StructureGraph:
    """Build a stable correlation graph via bootstrap selection."""
    if expr.ndim != 2:
        raise ValueError(f"expr must be 2D, got shape {expr.shape}.")
    n_samples, n_features = expr.shape
    if n_features <= 1 or k_max <= 0 or n_boot <= 0:
        return StructureGraph(
            neighbors=np.full((n_features, 0), -1, dtype=int),
            weights=np.zeros((n_features, 0), dtype=float),
        )
    if not (0.0 < subsample_frac <= 1.0):
        raise ValueError("subsample_frac must be in (0, 1].")

    rng = np.random.default_rng(seed)
    counts = [dict() for _ in range(n_features)]
    weight_sums = [dict() for _ in range(n_features)]

    subsample_n = max(2, int(round(subsample_frac * n_samples)))
    for _ in range(int(n_boot)):
        idx = rng.choice(n_samples, size=subsample_n, replace=False)
        graph = build_structure_graph(
            expr[idx],
            k=int(k_max),
            corr_method=corr_method,
            corr_estimator=corr_estimator,
            use_abs=use_abs,
            eps=eps,
        )
        for i in range(n_features):
            for j, w in zip(graph.neighbors[i], graph.weights[i]):
                if j < 0:
                    continue
                counts[i][int(j)] = counts[i].get(int(j), 0) + 1
                weight_sums[i][int(j)] = weight_sums[i].get(int(j), 0.0) + float(w)

    if min_freq <= 1.0:
        min_count = int(np.ceil(min_freq * n_boot))
    else:
        min_count = int(min_freq)
    min_count = max(min_count, 1)

    rows = []
    weights_rows = []
    for i in range(n_features):
        keep = [(j, counts[i][j], weight_sums[i][j]) for j in counts[i] if counts[i][j] >= min_count]
        keep.sort(key=lambda x: (x[1], abs(x[2])), reverse=True)
        if not keep:
            rows.append([])
            weights_rows.append([])
            continue
        keep = keep[: int(k_max)]
        nbrs = [j for j, _, _ in keep]
        wts = [weight_sums[i][j] / counts[i][j] for j, _, _ in keep]
        wts = np.asarray(wts, dtype=float)
        if use_abs:
            wts = np.abs(wts)
            denom = wts.sum()
        else:
            denom = np.sum(np.abs(wts))
        denom = denom if denom > eps else 1.0
        wts = wts / denom
        rows.append(nbrs)
        weights_rows.append(wts.tolist())

    k_out = max((len(r) for r in rows), default=0)
    if k_out == 0:
        return StructureGraph(
            neighbors=np.full((n_features, 0), -1, dtype=int),
            weights=np.zeros((n_features, 0), dtype=float),
        )

    neighbors = np.full((n_features, k_out), -1, dtype=int)
    weights = np.zeros((n_features, k_out), dtype=float)
    for i, (nbrs, wts) in enumerate(zip(rows, weights_rows)):
        if not nbrs:
            continue
        neighbors[i, : len(nbrs)] = np.asarray(nbrs, dtype=int)
        weights[i, : len(wts)] = np.asarray(wts, dtype=float)

    return StructureGraph(neighbors=neighbors, weights=weights)
