import numpy as np
import pandas as pd

from mass_sight import MassSightConfig
from mass_sight.multistudy import LoadedStudy, cluster_hub_mutual_top1


def _make_study(analysis_id: str, features: list[dict], expr: pd.DataFrame) -> LoadedStudy:
    feat = pd.DataFrame(features)
    feat["feature_id"] = feat["feature_id"].astype(str)
    feat = feat.set_index("feature_id", drop=False)
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    expr.columns = [str(c) for c in expr.columns]
    return LoadedStudy(analysis_id=analysis_id, features=feat, expr=expr)


def test_hub_clustering_full_join_keeps_unmatched_features_as_singletons():
    hub = _make_study(
        "HUB",
        [
            {"feature_id": "h1", "MZ": 100.0, "RT": 1.0, "Intensity": 10.0},
            {"feature_id": "h2", "MZ": 200.0, "RT": 2.0, "Intensity": 20.0},
        ],
        pd.DataFrame({"h1": [10.0, 0.0], "h2": [0.0, 20.0]}, index=["s1", "s2"]),
    )
    other = _make_study(
        "OTHER",
        [
            {"feature_id": "o1", "MZ": 100.0, "RT": 1.0, "Intensity": 11.0},  # matches h1
            {"feature_id": "o2", "MZ": 300.0, "RT": 3.0, "Intensity": 30.0},  # unmatched
        ],
        pd.DataFrame({"o1": [11.0, 0.0], "o2": [0.0, 30.0]}, index=["s1", "s2"]),
    )

    cfg = MassSightConfig(
        ppm=5.0,
        rt=0.0,
        tight_ppm=5.0,
        tight_rt=1.0,
        mz_shift_mode="none",
        allow_unmatched=True,
        rt_drift_model="none",
        ppm_drift_model="linear",
    )

    res = cluster_hub_mutual_top1([hub, other], cfg, hub_id="HUB", export_expr=True, include_unmatched=True)
    assert res.strategy == "hub_mutual_top1_full"

    singleton = "OTHER__o2"
    assert singleton in set(res.cluster_metadata["cluster_id"].astype(str))
    assert singleton in res.cluster_expr["HUB"].columns
    assert singleton in res.cluster_expr["OTHER"].columns

    # Structural missingness: singleton cluster doesn't exist in the hub study.
    assert res.cluster_expr["HUB"][singleton].isna().all()
    # Singleton cluster values come from the original unmatched feature.
    assert np.allclose(res.cluster_expr["OTHER"][singleton].to_numpy(dtype=float), other.expr["o2"].to_numpy(dtype=float))

    # Matched cluster: OTHER.o1 assigned into hub cluster h1.
    assert np.allclose(res.cluster_expr["OTHER"]["h1"].to_numpy(dtype=float), other.expr["o1"].to_numpy(dtype=float))
    # Hub-only cluster: h2 is structurally absent in OTHER.
    assert res.cluster_expr["OTHER"]["h2"].isna().all()
