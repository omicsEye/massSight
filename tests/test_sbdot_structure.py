import numpy as np
import pandas as pd

from mass_sight import MassSightConfig, match_features


def _count_correct(top1: pd.DataFrame, truth: dict[int, int]) -> int:
    correct = 0
    for _, row in top1.iterrows():
        i1 = int(row["id1"])
        i2 = int(row["id2"])
        if i2 < 0:
            continue
        if truth.get(i1) == i2:
            correct += 1
    return correct


def test_sbdot_structure_improves_on_ambiguous_local():
    # Local evidence is ambiguous and slightly favors a swapped match.
    study_a = pd.DataFrame(
        {
            "MZ": [100.0000, 100.0004, 150.0, 200.0],
            "RT": [1.0, 1.0, 2.0, 3.0],
            "Intensity": [1000.0, 900.0, 1200.0, 1100.0],
        }
    )
    study_b = pd.DataFrame(
        {
            "MZ": [100.0003, 100.0001, 150.0001, 200.0001],
            "RT": [1.0, 1.0, 2.0, 3.0],
            "Intensity": [1000.0, 900.0, 1200.0, 1100.0],
        }
    )

    # Correlation structure: (0,2) and (1,3) form two distinct clusters.
    x = np.linspace(1.0, 8.0, 8)
    y = np.array([1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0])
    expr_a = np.vstack([x, y, x + 0.01, y + 0.02]).T
    expr_b = np.vstack([x + 0.02, y + 0.01, x + 0.03, y + 0.03]).T

    truth = {0: 0, 1: 1, 2: 2, 3: 3}

    base_cfg = MassSightConfig(
        ppm=10.0,
        rt=0.1,
        tight_ppm=10.0,
        tight_rt=0.1,
        use_intensity=False,
        w_int=0.0,
        t_df=8.0,
        ot_epsilon=0.3,
        allow_unmatched=False,
    )
    base = match_features(study_a, study_b, base_cfg)
    base_correct = _count_correct(base.top1, truth)

    struct_cfg = MassSightConfig(
        ppm=10.0,
        rt=0.1,
        tight_ppm=10.0,
        tight_rt=0.1,
        use_intensity=False,
        w_int=0.0,
        t_df=8.0,
        ot_epsilon=0.3,
        allow_unmatched=False,
        use_structure=True,
        structure_alpha=5.0,
        structure_k=1,
        structure_iters=2,
        structure_min_samples=3,
    )
    struct = match_features(study_a, study_b, struct_cfg, expr_a=expr_a, expr_b=expr_b)
    struct_correct = _count_correct(struct.top1, truth)

    assert struct_correct > base_correct
