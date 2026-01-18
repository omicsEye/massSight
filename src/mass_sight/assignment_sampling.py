"""Sampling utilities for 1–1 matchings on sparse candidate graphs.

These functions are intended for *uncertainty/stability analysis* (e.g., perturb-and-match),
not for the primary MW benchmark evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class AssignmentSamplingConfig:
    n_samples: int = 100
    seed: int = 0
    gumbel_scale: float = 1.0
    null_utility: Optional[float] = None


@dataclass(frozen=True)
class ConsensusAlignmentConfig:
    """Consensus alignment from sampled 1–1 matchings.

    We sample multiple globally consistent 1–1 matchings (with explicit per-row nulls),
    compute edge inclusion frequencies ("support"), then return a single consensus
    1–1 assignment by solving a max-support matching. This yields a conservative
    alignment that keeps only edges that are stable under perturbations.

    `min_support` is a *stability* threshold (fraction of sampled matchings in which an
    edge appears), not a probability of correctness.
    """

    min_support: float = 0.25
    dummy_cost: float = 1e-6


def _gumbel(rng: np.random.Generator, size: int) -> np.ndarray:
    u = rng.random(size=size)
    u = np.clip(u, 1e-12, 1.0 - 1e-12)
    return -np.log(-np.log(u))


def sample_one_to_one_assignments(
    id1: np.ndarray,
    id2: np.ndarray,
    utility: np.ndarray,
    *,
    n1: int,
    n2: int,
    cfg: AssignmentSamplingConfig,
) -> np.ndarray:
    """Sample 1–1 matchings using Gumbel perturb-and-MAP assignment.

    The matching is sampled over a bipartite graph defined by (id1,id2) edges and their utilities.
    Unmatched rows are permitted by adding per-row dummy columns. The function returns an array
    `assignments` of shape (n_samples, n1) with values in {-1, 0..n2-1}.
    """
    id1 = np.asarray(id1, dtype=int)
    id2 = np.asarray(id2, dtype=int)
    utility = np.asarray(utility, dtype=float)
    if id1.shape != id2.shape or id1.shape != utility.shape:
        raise ValueError("id1/id2/utility must have the same shape.")
    if n1 <= 0 or n2 <= 0:
        raise ValueError("n1 and n2 must be positive.")

    n_samples = int(cfg.n_samples)
    if n_samples <= 0:
        raise ValueError("cfg.n_samples must be positive.")

    gumbel_scale = float(cfg.gumbel_scale)
    if not np.isfinite(gumbel_scale) or gumbel_scale <= 0:
        gumbel_scale = 1.0

    null_utility = cfg.null_utility
    if null_utility is None or not np.isfinite(float(null_utility)):
        finite = utility[np.isfinite(utility)]
        null_utility = float(np.nanquantile(finite, 0.2)) if finite.size else -10.0
    null_utility = float(null_utility)

    ok = (
        (id1 >= 0)
        & (id1 < int(n1))
        & (id2 >= 0)
        & (id2 < int(n2))
        & np.isfinite(utility)
    )
    id1 = id1[ok]
    id2 = id2[ok]
    utility = utility[ok]

    rng = np.random.default_rng(int(cfg.seed))
    assignments = np.full((n_samples, int(n1)), -1, dtype=int)

    from scipy.sparse import coo_matrix  # type: ignore
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching  # type: ignore

    dummy_cols = np.arange(int(n1), dtype=int) + int(n2)
    dummy_rows = np.arange(int(n1), dtype=int)

    for s in range(n_samples):
        noise_edges = gumbel_scale * _gumbel(rng, int(len(utility)))
        noise_dummy = gumbel_scale * _gumbel(rng, int(n1))
        cost_edges = -(utility + noise_edges)
        cost_dummy = -(null_utility + noise_dummy)

        rows = np.concatenate([id1, dummy_rows], axis=0)
        cols = np.concatenate([id2, dummy_cols], axis=0)
        data = np.concatenate([cost_edges, cost_dummy], axis=0)

        mat = coo_matrix((data, (rows, cols)), shape=(int(n1), int(n2) + int(n1))).tocsr()
        row_ind, col_ind = min_weight_full_bipartite_matching(mat)
        row_ind = np.asarray(row_ind, dtype=int)
        col_ind = np.asarray(col_ind, dtype=int)

        out = np.full(int(n1), -1, dtype=int)
        for r, c in zip(row_ind.tolist(), col_ind.tolist()):
            if c < int(n2):
                out[int(r)] = int(c)
        assignments[s, :] = out

    return assignments


def consensus_one_to_one_alignment(
    id1: np.ndarray,
    id2: np.ndarray,
    utility: np.ndarray,
    *,
    n1: int,
    n2: int,
    sample_cfg: AssignmentSamplingConfig,
    consensus_cfg: ConsensusAlignmentConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return a consensus 1–1 alignment and per-row support.

    Args:
        id1/id2/utility: candidate edges and their utilities (higher is better).
        n1/n2: number of features in dataset 1/2.
        sample_cfg: sampling config for Gumbel perturb-and-MAP.
        consensus_cfg: consensus config (support threshold + dummy cost).

    Returns:
        assignment: array shape (n1,) with values in {-1, 0..n2-1}
        support: array shape (n1,) giving the inclusion frequency of the chosen assignment;
                 for unmatched rows this is the frequency of being unmatched in the samples.
    """

    assigns = sample_one_to_one_assignments(
        id1,
        id2,
        utility,
        n1=int(n1),
        n2=int(n2),
        cfg=sample_cfg,
    )
    n_samples, n1_eff = assigns.shape
    if n1_eff != int(n1):
        raise ValueError("Internal error: sampled assignment shape mismatch.")

    min_support = float(consensus_cfg.min_support)
    if not np.isfinite(min_support):
        min_support = 0.7
    min_support = float(np.clip(min_support, 0.0, 1.0))

    dummy_cost = float(consensus_cfg.dummy_cost)
    if not np.isfinite(dummy_cost) or dummy_cost <= 0:
        dummy_cost = 1e-6

    # Compute edge inclusion frequencies from the sampled matchings.
    edge_counts: Dict[Tuple[int, int], int] = {}
    null_counts = np.zeros(int(n1), dtype=int)
    for s in range(int(n_samples)):
        row = assigns[s]
        for i in range(int(n1)):
            j = int(row[i])
            if j < 0:
                null_counts[i] += 1
            else:
                edge_counts[(i, j)] = int(edge_counts.get((i, j), 0)) + 1

    # Build a consensus matching on edges with support > min_support.
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    support_map: Dict[Tuple[int, int], float] = {}
    for (i, j), c in edge_counts.items():
        sup = float(c) / float(n_samples)
        support_map[(i, j)] = sup
        # Strictly require sup > min_support to avoid ties with the dummy option.
        if sup <= min_support + 1e-12:
            continue
        # Solve a min-cost matching, where cost is negative adjusted support.
        # adjusted_support = sup - min_support in (0,1].
        rows.append(int(i))
        cols.append(int(j))
        data.append(-(sup - min_support))

    # Always include per-row dummy columns (explicit null matches).
    dummy_rows = np.arange(int(n1), dtype=int)
    dummy_cols = int(n2) + np.arange(int(n1), dtype=int)
    rows.extend(dummy_rows.tolist())
    cols.extend(dummy_cols.tolist())
    data.extend([float(dummy_cost)] * int(n1))

    from scipy.sparse import coo_matrix  # type: ignore
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching  # type: ignore

    mat = coo_matrix((np.asarray(data, dtype=float), (np.asarray(rows, dtype=int), np.asarray(cols, dtype=int))), shape=(int(n1), int(n2) + int(n1))).tocsr()
    row_ind, col_ind = min_weight_full_bipartite_matching(mat)
    row_ind = np.asarray(row_ind, dtype=int)
    col_ind = np.asarray(col_ind, dtype=int)

    out = np.full(int(n1), -1, dtype=int)
    out_support = np.zeros(int(n1), dtype=float)
    for r, c in zip(row_ind.tolist(), col_ind.tolist()):
        if int(c) < int(n2):
            out[int(r)] = int(c)
            out_support[int(r)] = float(support_map.get((int(r), int(c)), 0.0))
        else:
            out[int(r)] = -1
            out_support[int(r)] = float(null_counts[int(r)]) / float(n_samples)

    return out, out_support


__all__ = [
    "AssignmentSamplingConfig",
    "sample_one_to_one_assignments",
    "ConsensusAlignmentConfig",
    "consensus_one_to_one_alignment",
]
