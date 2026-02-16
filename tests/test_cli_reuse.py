import json
from pathlib import Path

import pandas as pd

from mass_sight.cli import main as cli_main
from mass_sight.workbench import WorkbenchAnalysisBundle


def _write_mock_untarg(path: Path, *, mz_delta: float = 0.0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Samples\t100.0000_1.0000\t200.0000_2.0000",
        f"m/z\t{100.0000 + mz_delta:.4f}\t{200.0000 + mz_delta:.4f}",
        "rt\t1.0000\t2.0000",
        "S1\t1000\t2000",
        "S2\t1100\t1900",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_reuse_cli_groups_compatible_mw_analyses_and_writes_outputs(tmp_path, monkeypatch):
    def fake_fetch_workbench_analysis(
        analysis_id: str,
        *,
        cache_root: Path,
        timeout_s: float = 60.0,
        retries: int = 3,
        delay_s: float = 0.25,
        refresh: bool = False,
    ):
        aid = str(analysis_id).upper()
        mwtab_path = Path(cache_root) / "mwtab" / f"{aid}_mwtab.json"
        untarg_path = Path(cache_root) / "untarg" / f"{aid}_untarg.txt"
        mwtab_path.parent.mkdir(parents=True, exist_ok=True)
        mwtab_path.write_text("{}", encoding="utf-8")
        _write_mock_untarg(untarg_path, mz_delta=0.0001 if aid.endswith("2") else 0.0)
        metadata = {
            "analysis_id": aid,
            "study_id": "STX",
            "ion_mode_norm": "positive",
            "chromatography_norm": "reversed_phase",
            "instrument_type": "Orbitrap",
            "value_scale": "raw_like",
            "n_untarg_features": 2,
            "untarg_path": str(untarg_path),
            "mwtab_path": str(mwtab_path),
        }
        return WorkbenchAnalysisBundle(
            analysis_id=aid,
            mwtab_path=mwtab_path,
            untarg_path=untarg_path,
            metadata=metadata,
        )

    monkeypatch.setattr("mass_sight.cli.fetch_workbench_analysis", fake_fetch_workbench_analysis)

    out_dir = tmp_path / "reuse_out"
    rc = cli_main(
        [
            "reuse",
            "--analysis-id",
            "AN000001",
            "--analysis-id",
            "AN000002",
            "--out-dir",
            str(out_dir),
            "--strategy",
            "hub",
            "--use-intensity",
            "auto",
        ]
    )
    assert rc == 0

    metadata_tsv = out_dir / "analysis_metadata.tsv"
    groups_tsv = out_dir / "reuse_groups.tsv"
    summary_json = out_dir / "reuse_summary.json"
    run_manifest = out_dir / "run_manifest.json"
    group_dir = out_dir / "groups" / "positive__reversed_phase"
    cluster_map_tsv = group_dir / "cluster_map.tsv"

    assert metadata_tsv.exists()
    assert groups_tsv.exists()
    assert summary_json.exists()
    assert run_manifest.exists()
    assert cluster_map_tsv.exists()

    metadata_df = pd.read_csv(metadata_tsv, sep="\t")
    assert set(metadata_df["analysis_id"].astype(str)) == {"AN000001", "AN000002"}

    groups_df = pd.read_csv(groups_tsv, sep="\t")
    assert int(groups_df.loc[0, "n_studies"]) == 2
    assert str(groups_df.loc[0, "group_key"]) == "positive__reversed_phase"
    assert str(groups_df.loc[0, "use_intensity_requested"]) == "auto"

    payload = json.loads(summary_json.read_text(encoding="utf-8"))
    assert int(payload.get("n_requested_analyses", 0)) == 2
    assert int(payload.get("n_fetched_analyses", 0)) == 2
    assert int(payload.get("n_groups", 0)) == 1
