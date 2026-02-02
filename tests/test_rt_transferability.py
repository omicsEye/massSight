import numpy as np
import pandas as pd

from mass_sight import MassSightConfig, match_features


def test_rt_unit_normalization_seconds_to_minutes_is_recorded():
    n = 40
    mz = np.linspace(100.0, 500.0, n)
    rt_min = np.linspace(1.0, 20.0, n)
    rt_sec = rt_min * 60.0
    study_a = pd.DataFrame({"MZ": mz, "RT": rt_min, "Intensity": np.full(n, 1e6)})
    study_b = pd.DataFrame({"MZ": mz, "RT": rt_sec, "Intensity": np.full(n, 1e6)})

    cfg = MassSightConfig(ppm=5.0, rt=0.0, use_intensity=False, w_int=0.0, allow_unmatched=False)
    res = match_features(study_a, study_b, cfg)
    assert abs(float(res.drift_params.get("rt_unit_scale_study_b", 1.0)) - (1.0 / 60.0)) < 1e-12


def test_rt_transferability_fits_linear_scale_offset_and_marks_policy_linear():
    n = 60
    mz = np.linspace(100.0, 900.0, n)
    rt_a = np.linspace(0.5, 25.0, n)
    # Dataset B is a linear transform of A (minutes) plus tiny noise.
    rt_b = 2.0 * rt_a + 1.0 + np.linspace(-0.01, 0.01, n)
    study_a = pd.DataFrame({"MZ": mz, "RT": rt_a, "Intensity": np.full(n, 1e6)})
    study_b = pd.DataFrame({"MZ": mz, "RT": rt_b, "Intensity": np.full(n, 1e6)})

    cfg = MassSightConfig(
        ppm=5.0,
        rt=0.0,
        use_intensity=False,
        w_int=0.0,
        allow_unmatched=False,
        rt_transfer_anchor_ppm=2.0,
        rt_transfer_min_anchors=25,
    )
    res = match_features(study_a, study_b, cfg)
    assert res.drift_params.get("rt_policy") == "linear"
    assert abs(float(res.drift_params.get("rt_scale", 1.0)) - 0.5) < 1e-2
    assert abs(float(res.drift_params.get("rt_offset", 0.0)) - (-0.5)) < 1e-1

