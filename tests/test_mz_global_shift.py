import numpy as np
import pandas as pd

from mass_sight import MassSightConfig, match_features
from mass_sight.lcms_utils import PROTON_MASS


def test_mz_global_shift_applies_proton_shift_for_neutral_vs_ion_positive():
    n = 80
    mz_neutral = np.linspace(100.0, 900.0, n)
    rt = np.linspace(1.0, 20.0, n)

    # Study A reports neutral masses; study B reports [M+H]+ ion m/z.
    study_a = pd.DataFrame({"MZ": mz_neutral, "RT": rt, "Intensity": np.full(n, 1e6)})
    study_b = pd.DataFrame({"MZ": mz_neutral + float(PROTON_MASS), "RT": rt, "Intensity": np.full(n, 1e6)})

    cfg = MassSightConfig(
        ppm=5.0,
        rt=0.0,
        tight_ppm=2.0,
        tight_rt=0.0,
        allow_unmatched=False,
        polarity="positive",
        mz_global_shift_model="auto",
        mz_semantics_study_a="neutral_mass",
        mz_semantics_study_b="ion_mz",
        mz_global_shift_anchor_ppm=2.0,
        mz_global_shift_min_anchors=10,
        mz_global_shift_min_gain=2,
    )

    res = match_features(study_a, study_b, cfg)
    assert bool(res.drift_params.get("mz_shift_applied", False)) is True

    # The alignment may choose to shift either A (+H) or B (-H); accept either.
    sa = float(res.drift_params.get("mz_shift_da_study_a", 0.0) or 0.0)
    sb = float(res.drift_params.get("mz_shift_da_study_b", 0.0) or 0.0)
    ok = (abs(sa - float(PROTON_MASS)) < 1e-6 and abs(sb) < 1e-12) or (
        abs(sb + float(PROTON_MASS)) < 1e-6 and abs(sa) < 1e-12
    )
    assert ok

    top = res.top1
    assert int(len(top)) == n
    assert (top["id2"].to_numpy(dtype=int) >= 0).all()
    assert np.all(top["id2"].to_numpy(dtype=int) == np.arange(n, dtype=int))

