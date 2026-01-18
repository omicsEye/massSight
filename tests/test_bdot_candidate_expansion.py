import numpy as np
import pandas as pd

from mass_sight import MassSightConfig, match_features
from mass_sight.lcms_utils import PROTON_MASS


def _toy_shifted_pair(n: int = 8, *, delta: float = PROTON_MASS) -> tuple[pd.DataFrame, pd.DataFrame]:
    mz1 = np.linspace(100.0, 400.0, n)
    rt1 = np.linspace(1.0, 10.0, n)
    study_a = pd.DataFrame(
        {
            "MZ": mz1,
            "RT": rt1,
            "Intensity": np.full(n, 1e6),
        }
    )
    study_b = pd.DataFrame(
        {
            "MZ": mz1 + float(delta),
            "RT": rt1,
            "Intensity": np.full(n, 1e6),
        }
    )
    return study_a, study_b


def test_bdot_candidate_expansion_none_blocks_shifted_matches():
    study_a, study_b = _toy_shifted_pair(n=6, delta=PROTON_MASS)
    cfg = MassSightConfig(
        ppm=10.0,
        rt=0.0,
        use_intensity=False,
        w_int=0.0,
        mz_shift_mode="none",
    )
    res = match_features(study_a, study_b, cfg)
    assert res.candidates.empty


def test_bdot_candidate_expansion_manual_recovers_candidates():
    study_a, study_b = _toy_shifted_pair(n=6, delta=PROTON_MASS)
    cfg = MassSightConfig(
        ppm=10.0,
        rt=0.0,
        tight_ppm=10.0,
        tight_rt=0.5,
        use_intensity=False,
        w_int=0.0,
        mz_shift_mode="manual",
        mz_shift_manual_deltas_da=[PROTON_MASS],
        mz_shift_penalty=0.0,
    )
    res = match_features(study_a, study_b, cfg)
    assert not res.candidates.empty
    assert "mz_shift_da" in res.candidates.columns
    assert np.allclose(res.candidates["mz_shift_da"].to_numpy(dtype=float), PROTON_MASS, atol=1e-9)

    # Shift-aware ppm residual should be near zero for true pairs.
    assert np.nanmax(np.abs(res.candidates["delta_ppm_raw"].to_numpy(dtype=float))) < 1e-6

    top1 = res.top1.sort_values("id1").reset_index(drop=True)
    assert np.array_equal(top1["id1"].to_numpy(dtype=int), np.arange(len(study_a)))
    assert np.array_equal(top1["id2"].to_numpy(dtype=int), np.arange(len(study_a)))


def test_bdot_candidate_expansion_dedup_picks_best_shift():
    study_a, study_b = _toy_shifted_pair(n=1, delta=PROTON_MASS)
    cfg = MassSightConfig(
        ppm=150.0,  # wide enough so both 1.0 and PROTON_MASS neighborhoods include the point
        rt=0.0,
        tight_ppm=150.0,
        tight_rt=1.0,
        use_intensity=False,
        w_int=0.0,
        mz_shift_mode="manual",
        mz_shift_manual_deltas_da=[1.0, PROTON_MASS],
        mz_shift_penalty=0.0,
    )
    res = match_features(study_a, study_b, cfg)
    assert len(res.candidates) == 1
    assert abs(float(res.candidates.iloc[0]["mz_shift_da"]) - float(PROTON_MASS)) < 1e-9


def test_bdot_candidate_expansion_auto_selects_supported_shift():
    study_a, study_b = _toy_shifted_pair(n=30, delta=PROTON_MASS)
    cfg = MassSightConfig(
        ppm=10.0,
        rt=0.0,
        tight_ppm=10.0,
        tight_rt=0.5,
        use_intensity=False,
        w_int=0.0,
        mz_shift_mode="auto",
        mz_shift_top_k=3,
        mz_shift_min_support_count=25,
        mz_shift_min_support_frac=0.5,
        mz_shift_penalty=0.0,
    )
    res = match_features(study_a, study_b, cfg)
    sel = res.drift_params.get("mz_shift_selected_deltas_da")
    assert isinstance(sel, list)
    assert any(abs(float(d) - float(PROTON_MASS)) < 1e-6 for d in sel)
