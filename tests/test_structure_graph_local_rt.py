import numpy as np

from mass_sight.structure_graph import build_structure_graph_local_rt


def test_structure_graph_local_rt_respects_rt_window():
    # Two RT-local correlated pairs: (0,1) and (2,3). Feature 4 is isolated.
    rt = np.array([0.00, 0.01, 0.50, 0.51, 1.00], dtype=float)
    x = np.linspace(0.0, 1.0, 20)
    y = np.sin(np.linspace(0.0, 3.0, 20))
    expr = np.vstack(
        [
            x,
            x + 1e-3,
            y,
            y + 1e-3,
            np.random.default_rng(0).normal(size=20),
        ]
    ).T  # (n_samples, n_features)

    g = build_structure_graph_local_rt(expr, rt, k=1, rt_window=0.05, corr_method="pearson", use_abs=True)

    # Feature 0 should see only feature 1 within the RT window.
    assert int(g.neighbors[0, 0]) == 1
    # Feature 2 should see only feature 3 within the RT window.
    assert int(g.neighbors[2, 0]) == 3
    # Feature 4 has no neighbor within 0.05 minutes.
    assert int(g.neighbors[4, 0]) == -1

