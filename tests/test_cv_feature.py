import pandas as pd

from mass_sight import MassSightConfig, match_features


def test_cv_evidence_prefers_similar_cv():
    study_a = pd.DataFrame(
        {
            "MZ": [100.0],
            "RT": [1.0],
            "Intensity": [1000.0],
            "CV": [0.10],
        }
    )
    study_b = pd.DataFrame(
        {
            "MZ": [100.0002, 100.0002],
            "RT": [1.0, 1.0],
            "Intensity": [1000.0, 1000.0],
            "CV": [0.10, 0.50],
        }
    )

    cfg = MassSightConfig(
        ppm=50.0,
        rt=0.0,
        tight_ppm=10.0,
        tight_rt=0.0,
        allow_unmatched=False,
        w_ppm=0.0,
        w_rt=0.0,
        w_roi=0.0,
        use_intensity=False,
        w_int=0.0,
        use_cv=True,
        w_cv=2.0,
        ot_mode="semi_relaxed",
    )
    res = match_features(study_a, study_b, cfg)
    top = res.top1.sort_values("id1").reset_index(drop=True)
    assert int(top.loc[0, "id2"]) == 0

    cand = res.candidates
    assert not cand.empty
    s = cand.loc[cand["id1"] == 0].sort_values("id2").reset_index(drop=True)
    assert float(s.loc[0, "loglik_local"]) > float(s.loc[1, "loglik_local"])

