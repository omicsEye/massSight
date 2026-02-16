import json

import pandas as pd

from mass_sight.cli import main as cli_main


def test_match_cli_writes_run_manifest_by_default(tmp_path):
    study_a = pd.DataFrame(
        {
            "MZ": [100.0, 200.0, 300.0],
            "RT": [1.0, 2.0, 3.0],
            "Intensity": [1000.0, 2000.0, 3000.0],
        }
    )
    study_b = pd.DataFrame(
        {
            "MZ": [100.0001, 200.0001, 300.0001],
            "RT": [1.0, 2.0, 3.0],
            "Intensity": [1100.0, 2100.0, 3100.0],
        }
    )

    a_path = tmp_path / "study_a.csv"
    b_path = tmp_path / "study_b.csv"
    out_candidates = tmp_path / "candidates.csv"
    out_top1 = tmp_path / "top1.csv"
    out_manifest = tmp_path / "candidates.run_manifest.json"

    study_a.to_csv(a_path, index=False)
    study_b.to_csv(b_path, index=False)

    rc = cli_main(
        [
            "match",
            str(a_path),
            str(b_path),
            "--out-candidates",
            str(out_candidates),
            "--out-top1",
            str(out_top1),
        ]
    )

    assert rc == 0
    assert out_candidates.exists()
    assert out_top1.exists()
    assert out_manifest.exists()

    payload = json.loads(out_manifest.read_text())
    assert payload.get("command") == "match"
    assert payload.get("mass_sight_version")
    assert int(payload.get("n_features_study_a", 0)) == 3
    assert int(payload.get("n_features_study_b", 0)) == 3
    assert int(payload.get("n_candidates", -1)) >= 1
    assert int(payload.get("n_top1_rows", -1)) == 3
