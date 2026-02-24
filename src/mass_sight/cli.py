import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .matcher import MassSightConfig, match_features
from .multistudy import (
    cluster_hub_barycenter_ot,
    cluster_hub_mutual_top1,
    cluster_symmetric_mutual_graph,
    load_manifest,
    load_study,
)
from .workbench import (
    bundle_to_study_spec,
    fetch_workbench_analysis,
    fetch_workbench_study_data,
    group_workbench_analyses,
    extract_mwtab_labeled_metabolites,
    normalize_metabolite_name,
    parse_workbench_study_data_json,
    resolve_use_intensity_mode,
)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _default_manifest_path(args: argparse.Namespace) -> Path:
    if getattr(args, "out_manifest", None):
        return Path(str(args.out_manifest))
    if str(args.cmd) == "match":
        out_candidates = Path(str(args.out_candidates))
        return out_candidates.parent / f"{out_candidates.stem}.run_manifest.json"
    if str(args.cmd) in {"cluster", "reuse"}:
        return Path(str(args.out_dir)) / "run_manifest.json"
    return Path("run_manifest.json")


def _serialize_args(args: argparse.Namespace) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in vars(args).items():
        if isinstance(value, Path):
            out[key] = str(value)
        else:
            out[key] = value
    return out


def _write_run_manifest(
    *,
    args: argparse.Namespace,
    started_at: str,
    ended_at: str,
    runtime_seconds: float,
    outputs: list[str],
    extras: dict[str, Any],
) -> Path:
    out_path = _default_manifest_path(args)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "command": str(args.cmd),
        "mass_sight_version": str(__version__),
        "started_at_utc": str(started_at),
        "ended_at_utc": str(ended_at),
        "runtime_seconds": float(runtime_seconds),
        "python_version": str(sys.version.split()[0]),
        "platform": str(platform.platform()),
        "parameters": _serialize_args(args),
        "outputs": [str(p) for p in outputs],
        **extras,
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return out_path


def _add_match_parser(sub):
    p = sub.add_parser("match", help="Match two studies with massSight")
    p.add_argument("study_a", type=str, help="Path to study_a feature table CSV")
    p.add_argument("study_b", type=str, help="Path to study_b feature table CSV")
    p.add_argument("--out-candidates", required=True, type=str, help="Output CSV for candidate edges")
    p.add_argument("--out-top1", default=None, type=str, help="Optional output CSV for top-1 matches")

    # Core windows
    p.add_argument("--ppm", type=float, default=20.0, help="m/z tolerance in ppm")
    p.add_argument("--rt", type=float, default=0.0, help="RT candidate window in minutes (0 disables RT gating)")
    p.add_argument("--tight-ppm", type=float, default=7.0, help="Tight ppm window for drift fitting")
    p.add_argument("--tight-rt", type=float, default=0.0, help="Tight RT window for drift fitting")

    # Intensity / unmatched mass
    p.add_argument("--use-intensity", dest="use_intensity", action="store_true", help="Enable intensity residual term")
    p.add_argument("--no-intensity", dest="use_intensity", action="store_false",
                   help="Disable intensity residual term (default)")
    p.set_defaults(use_intensity=False)
    p.add_argument(
        "--allow-unmatched",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add a null column to allow 'no match' (default: enabled)",
    )
    p.add_argument("--null-mass", type=float, default=0.1, help="OT mass assigned to the null column when enabled")

    # Schema
    p.add_argument("--mz-col", type=str, default="MZ")
    p.add_argument("--rt-col", type=str, default="RT")
    p.add_argument("--intensity-col", type=str, default="Intensity")
    p.add_argument("--polarity", choices=["positive", "negative"], default="positive")
    p.add_argument(
        "--mz-semantics-a",
        choices=["auto", "ion_mz", "neutral_mass"],
        default="auto",
        help="Interpretation of study_a MZ column (default: auto).",
    )
    p.add_argument(
        "--mz-semantics-b",
        choices=["auto", "ion_mz", "neutral_mass"],
        default="auto",
        help="Interpretation of study_b MZ column (default: auto).",
    )
    p.add_argument(
        "--mz-global-shift-model",
        choices=["auto", "none", "proton_only", "common_adducts"],
        default="auto",
        help="Global mass-axis alignment to handle neutral-mass feature tables (default: auto).",
    )
    p.add_argument("--mz-global-shift-anchor-ppm", type=float, default=5.0, help="Anchor ppm for m/z shift inference.")
    p.add_argument("--mz-global-shift-min-anchors", type=int, default=25, help="Min MNN anchors to apply a global shift.")
    p.add_argument("--mz-global-shift-min-gain", type=int, default=5, help="Min anchor-count gain over 0-shift to apply.")
    p.add_argument("--mz-global-shift-max-eval", type=int, default=5000, help="Max m/z values to use per study for shift inference.")
    p.add_argument(
        "--rt-transfer-min-inlier-frac",
        type=float,
        default=0.02,
        help="RT transferability: minimum inlier fraction for declaring RT transferable (default: 0.02).",
    )
    p.add_argument(
        "--retention-nontransferable-policy",
        choices=["disable", "soft"],
        default="soft",
        help="When RT is non-transferable, either disable retention or use it as a weak tie-breaker (default: soft).",
    )
    p.add_argument(
        "--retention-nontransferable-weight",
        type=float,
        default=0.25,
        help="Base weight multiplier for retention under non-transferable RT when policy=soft (default: 0.25).",
    )
    p.add_argument(
        "--capture-manifest",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write run metadata manifest with version, parameters, runtime, and output paths (default: enabled).",
    )
    p.add_argument(
        "--out-manifest",
        type=str,
        default=None,
        help="Optional path for run manifest JSON (default: <out-candidates-stem>.run_manifest.json).",
    )
    return p


def _add_cluster_parser(sub):
    p = sub.add_parser("cluster", help="Cluster multiple untargeted studies with massSight")
    p.add_argument("--manifest", required=True, type=str, help="Path to YAML/JSON multi-study manifest")
    p.add_argument("--out-dir", required=True, type=str, help="Output directory")

    p.add_argument(
        "--strategy",
        type=str,
        default="both",
        choices=["hub", "symmetric", "both"],
        help="Clustering strategy: hub (primary), symmetric (sensitivity), or both (default).",
    )
    p.add_argument("--hub-id", type=str, default=None, help="Optional hub study ID for hub strategy")
    p.add_argument(
        "--hub-alignment",
        choices=["auto", "mutual_top1", "barycenter"],
        default="auto",
        help="Hub strategy alignment: auto (barycenter for k>=3), mutual_top1, or barycenter.",
    )
    p.add_argument(
        "--full-join",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Include unmatched features from all studies as singleton clusters (full outer join).",
    )
    p.add_argument("--no-export-expr", dest="export_expr", action="store_false", default=True,
                   help="Do not export sample×cluster matrices (default: export when available).")

    # massSight config (minimal set; advanced tuning via Python API)
    p.add_argument("--ppm", type=float, default=20.0, help="Initial ppm window for candidate search")
    p.add_argument("--rt", type=float, default=0.0, help="RT gate during candidate generation (0 disables)")
    p.add_argument("--tight-ppm", type=float, default=7.0, help="Tight ppm window for drift fitting")
    p.add_argument("--tight-rt", type=float, default=0.0, help="Tight RT window for drift fitting")
    p.add_argument("--polarity", choices=["positive", "negative"], default=None,
                   help="Override polarity (default: infer/validate from manifest).")
    p.add_argument(
        "--mz-global-shift-model",
        choices=["auto", "none", "proton_only", "common_adducts"],
        default="auto",
        help="Global mass-axis alignment to handle neutral-mass feature tables (default: auto).",
    )
    p.add_argument("--mz-global-shift-anchor-ppm", type=float, default=5.0, help="Anchor ppm for m/z shift inference.")
    p.add_argument("--mz-global-shift-min-anchors", type=int, default=25, help="Min MNN anchors to apply a global shift.")
    p.add_argument("--mz-global-shift-min-gain", type=int, default=5, help="Min anchor-count gain over 0-shift to apply.")
    p.add_argument("--mz-global-shift-max-eval", type=int, default=5000, help="Max m/z values to use per study for shift inference.")

    p.add_argument(
        "--drift-fit-bootstrap",
        action="store_true",
        default=False,
        help="Fit drift models on a bootstrap resample of tight matches (uncertainty propagation).",
    )
    p.add_argument("--drift-bootstrap-seed", type=int, default=0, help="RNG seed for drift bootstrap resampling.")
    p.add_argument(
        "--capture-manifest",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write run metadata manifest with version, parameters, runtime, and output paths (default: enabled).",
    )
    p.add_argument(
        "--out-manifest",
        type=str,
        default=None,
        help="Optional path for run manifest JSON (default: <out-dir>/run_manifest.json).",
    )
    return p


def _add_reuse_parser(sub):
    p = sub.add_parser(
        "reuse",
        help="Fetch Metabolomics Workbench analyses and run grouped cross-study reuse clustering",
    )
    p.add_argument(
        "--analysis-id",
        dest="analysis_ids",
        action="append",
        default=[],
        help="Workbench analysis ID (repeat flag for multiple analyses, e.g. AN000001).",
    )
    p.add_argument(
        "--analysis-file",
        type=str,
        default=None,
        help="Optional text file with one analysis ID per line, or a JSON selection manifest from `mass_sight find`.",
    )
    p.add_argument("--out-dir", required=True, type=str, help="Output directory")
    p.add_argument(
        "--cache-dir",
        type=str,
        default=None,
        help="Optional cache directory for fetched mwTab and untarg files (default: <out-dir>/cache).",
    )
    p.add_argument(
        "--refresh-fetch",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Re-fetch Workbench files even when cached locally (default: disabled).",
    )
    p.add_argument("--fetch-timeout-s", type=float, default=60.0, help="HTTP timeout in seconds for MW fetches.")
    p.add_argument("--fetch-retries", type=int, default=3, help="Retries per MW endpoint request.")
    p.add_argument("--fetch-delay-s", type=float, default=0.25, help="Delay between MW requests in seconds.")

    p.add_argument("--min-group-size", type=int, default=2, help="Minimum analyses per compatibility group.")
    p.add_argument(
        "--allow-unknown-strata",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Allow groups with unknown polarity/chromatography (default: disabled).",
    )
    p.add_argument(
        "--export-labeled-metabolites",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Export per-analysis labeled metabolite tables and shared-labeled-metabolite overlaps (default: enabled).",
    )
    p.add_argument(
        "--fetch-targeted-data",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Fetch Workbench study-level named-metabolite matrices (/study_id/.../data) for targeted/meta-analysis (default: disabled).",
    )
    p.add_argument(
        "--targeted-max-mib",
        type=float,
        default=50.0,
        help="Max MiB to download per study_id targeted data payload when --fetch-targeted-data is enabled.",
    )
    p.add_argument(
        "--use-intensity",
        dest="use_intensity_mode",
        choices=["off", "auto", "on"],
        default="off",
        help=(
            "Control use of intensity similarity across studies: "
            "'off' (ignore intensity), 'auto' (enable only when metadata suggests comparable raw-like values), "
            "'on' (always include intensity)."
        ),
    )

    p.add_argument(
        "--strategy",
        type=str,
        default="both",
        choices=["hub", "symmetric", "both"],
        help="Clustering strategy per compatibility group.",
    )
    p.add_argument(
        "--hub-alignment",
        choices=["auto", "mutual_top1", "barycenter"],
        default="auto",
        help="Hub alignment for hub strategy: auto, mutual_top1, or barycenter.",
    )
    p.add_argument(
        "--full-join",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Include unmatched features as singleton clusters (full outer join).",
    )
    p.add_argument("--no-export-expr", dest="export_expr", action="store_false", default=True)

    p.add_argument("--ppm", type=float, default=20.0, help="Initial ppm window for candidate search")
    p.add_argument("--rt", type=float, default=0.0, help="RT gate during candidate generation (0 disables)")
    p.add_argument("--tight-ppm", type=float, default=7.0, help="Tight ppm window for drift fitting")
    p.add_argument("--tight-rt", type=float, default=0.0, help="Tight RT window for drift fitting")
    p.add_argument(
        "--polarity",
        choices=["positive", "negative"],
        default=None,
        help="Force polarity for all groups (otherwise inferred from Workbench ion mode metadata).",
    )
    p.add_argument(
        "--mz-global-shift-model",
        choices=["auto", "none", "proton_only", "common_adducts"],
        default="auto",
        help="Global mass-axis alignment to handle neutral-mass feature tables (default: auto).",
    )
    p.add_argument("--mz-global-shift-anchor-ppm", type=float, default=5.0, help="Anchor ppm for m/z shift inference.")
    p.add_argument("--mz-global-shift-min-anchors", type=int, default=25, help="Min MNN anchors to apply a global shift.")
    p.add_argument("--mz-global-shift-min-gain", type=int, default=5, help="Min anchor-count gain over 0-shift to apply.")
    p.add_argument("--mz-global-shift-max-eval", type=int, default=5000, help="Max m/z values to use per study for shift inference.")

    p.add_argument(
        "--drift-fit-bootstrap",
        action="store_true",
        default=False,
        help="Fit drift models on a bootstrap resample of tight matches.",
    )
    p.add_argument("--drift-bootstrap-seed", type=int, default=0, help="RNG seed for drift bootstrap resampling.")

    p.add_argument(
        "--capture-manifest",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write run metadata manifest with version, parameters, runtime, and output paths (default: enabled).",
    )
    p.add_argument(
        "--out-manifest",
        type=str,
        default=None,
        help="Optional path for run manifest JSON (default: <out-dir>/run_manifest.json).",
    )
    return p


def _add_find_parser(sub):
    p = sub.add_parser("find", help="Interactive TUI to find Workbench studies/analyses (start from disease)")
    p.add_argument("--out", required=True, type=str, help="Output selection manifest JSON path.")
    p.add_argument(
        "--cache-dir",
        type=str,
        default=None,
        help="Optional cache directory for fetched Workbench indices (default: <out-dir>/cache).",
    )
    p.add_argument(
        "--refresh-fetch",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Re-fetch Workbench indices even when cached locally (default: disabled).",
    )
    p.add_argument("--fetch-timeout-s", type=float, default=60.0, help="HTTP timeout in seconds for MW fetches.")
    p.add_argument("--fetch-retries", type=int, default=3, help="Retries per MW endpoint request.")
    p.add_argument("--fetch-delay-s", type=float, default=0.25, help="Delay between MW requests in seconds.")
    return p


def _collect_analysis_ids(args: argparse.Namespace) -> list[str]:
    values: list[str] = []
    values.extend([str(x).strip() for x in (args.analysis_ids or []) if str(x).strip()])
    if getattr(args, "analysis_file", None):
        path = Path(str(args.analysis_file))
        if not path.exists():
            raise FileNotFoundError(path)
        raw_text = path.read_text(encoding="utf-8")

        # Support selection manifests from `mass_sight find` in addition to plain text lists.
        # If parsing fails, fall back to line-based extraction.
        parsed_json = None
        try:
            parsed_json = json.loads(raw_text)
        except Exception:
            parsed_json = None

        if isinstance(parsed_json, dict) and isinstance(parsed_json.get("analysis_ids"), list):
            for item in parsed_json.get("analysis_ids") or []:
                s = str(item).strip()
                if s:
                    values.append(s)
        elif isinstance(parsed_json, list):
            for item in parsed_json:
                s = str(item).strip()
                if s:
                    values.append(s)
        else:
            for line in raw_text.splitlines():
                raw = str(line).strip()
                if not raw or raw.startswith("#"):
                    continue
                values.append(raw)
    deduped: list[str] = []
    seen: set[str] = set()
    for raw in values:
        s = str(raw).upper().strip()
        if not s:
            continue
        if not s.startswith("AN"):
            try:
                s = f"AN{int(s):06d}"
            except Exception:
                pass
        if s in seen:
            continue
        seen.add(s)
        deduped.append(s)
    return deduped


def _run_cluster_strategies(
    *,
    studies,
    cfg: MassSightConfig,
    out_dir: Path,
    strategy: str,
    hub_id: str | None,
    hub_alignment: str,
    full_join: bool,
    export_expr: bool,
) -> tuple[list[str], list[str]]:
    strategies = [strategy] if strategy != "both" else ["hub", "symmetric"]
    outputs: list[str] = []
    for strat in strategies:
        if strat == "hub":
            hub_alignment_norm = str(hub_alignment).lower().strip()
            use_barycenter = hub_alignment_norm == "barycenter" or (
                hub_alignment_norm == "auto" and len(studies) >= 3
            )
            if use_barycenter:
                res = cluster_hub_barycenter_ot(
                    studies,
                    cfg,
                    hub_id=hub_id,
                    include_unmatched=bool(full_join),
                    export_expr=bool(export_expr),
                )
            else:
                res = cluster_hub_mutual_top1(
                    studies,
                    cfg,
                    hub_id=hub_id,
                    export_expr=bool(export_expr),
                    include_unmatched=bool(full_join),
                )
        else:
            res = cluster_symmetric_mutual_graph(studies, cfg, export_expr=bool(export_expr))

        sub_out = out_dir / res.strategy if strategy == "both" else out_dir
        sub_out.mkdir(parents=True, exist_ok=True)

        cluster_metadata_path = sub_out / "cluster_metadata.tsv"
        cluster_map_path = sub_out / "cluster_map.tsv"
        diagnostics_path = sub_out / "diagnostics.json"
        res.cluster_metadata.to_csv(cluster_metadata_path, sep="\t", index=False)
        res.cluster_map.to_csv(cluster_map_path, sep="\t", index=False)
        outputs.extend([str(cluster_metadata_path), str(cluster_map_path), str(diagnostics_path)])
        for study_id, mat in res.cluster_expr.items():
            expr_path = sub_out / f"cluster_expr_{study_id}.tsv.gz"
            mat.to_csv(expr_path, sep="\t", compression="gzip")
            outputs.append(str(expr_path))
        diagnostics_path.write_text(json.dumps(res.diagnostics, indent=2, sort_keys=True))

    return outputs, [str(x) for x in strategies]


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    ap = argparse.ArgumentParser(prog="mass_sight", description="massSight matching for LC–MS")
    sub = ap.add_subparsers(dest="cmd", required=True)
    _add_match_parser(sub)
    _add_cluster_parser(sub)
    _add_reuse_parser(sub)
    _add_find_parser(sub)
    args = ap.parse_args(argv)
    started_at = _utc_now_iso()
    t0 = time.perf_counter()

    if args.cmd == "match":
        study_a = pd.read_csv(args.study_a)
        study_b = pd.read_csv(args.study_b)
        cfg = MassSightConfig(
            ppm=float(args.ppm),
            rt=float(args.rt),
            tight_ppm=float(args.tight_ppm),
            tight_rt=float(args.tight_rt),
            use_intensity=bool(args.use_intensity),
            allow_unmatched=bool(args.allow_unmatched),
            null_mass=float(args.null_mass),
            mz_col=str(args.mz_col),
            rt_col=str(args.rt_col),
            intensity_col=str(args.intensity_col) if args.intensity_col else None,
            polarity=str(args.polarity),
            mz_semantics_study_a=str(getattr(args, "mz_semantics_a", "auto")),
            mz_semantics_study_b=str(getattr(args, "mz_semantics_b", "auto")),
            mz_global_shift_model=str(getattr(args, "mz_global_shift_model", "auto")),
            mz_global_shift_anchor_ppm=float(getattr(args, "mz_global_shift_anchor_ppm", 5.0)),
            mz_global_shift_min_anchors=int(getattr(args, "mz_global_shift_min_anchors", 25)),
            mz_global_shift_min_gain=int(getattr(args, "mz_global_shift_min_gain", 5)),
            mz_global_shift_max_eval=int(getattr(args, "mz_global_shift_max_eval", 5000)),
            rt_transfer_min_inlier_frac=float(getattr(args, "rt_transfer_min_inlier_frac", 0.02)),
            retention_nontransferable_policy=str(getattr(args, "retention_nontransferable_policy", "soft")),
            retention_nontransferable_weight=float(getattr(args, "retention_nontransferable_weight", 0.25)),
        )
        res = match_features(study_a, study_b, cfg)
        res.candidates.to_csv(args.out_candidates, index=False)
        outputs = [str(args.out_candidates)]
        if args.out_top1:
            res.top1.to_csv(args.out_top1, index=False)
            outputs.append(str(args.out_top1))
        if bool(args.capture_manifest):
            ended_at = _utc_now_iso()
            runtime_seconds = time.perf_counter() - t0
            manifest_path = _write_run_manifest(
                args=args,
                started_at=started_at,
                ended_at=ended_at,
                runtime_seconds=runtime_seconds,
                outputs=outputs,
                extras={
                    "n_features_study_a": int(len(study_a)),
                    "n_features_study_b": int(len(study_b)),
                    "n_candidates": int(len(res.candidates)),
                    "n_top1_rows": int(len(res.top1)),
                },
            )
            outputs.append(str(manifest_path))
        print(
            f"wrote {args.out_candidates} ({len(res.candidates)} candidates); "
            f"top1 rows={len(res.top1)}"
        )
        return 0

    if args.cmd == "cluster":
        manifest_path = Path(args.manifest)
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        specs = load_manifest(manifest_path)
        studies = [load_study(s) for s in specs]

        manifest_polarities = sorted({str(s.spec.polarity).lower().strip() for s in studies if s.spec is not None})
        if args.polarity:
            polarity = str(args.polarity).lower().strip()
        elif len(manifest_polarities) == 1:
            polarity = manifest_polarities[0]
        else:
            raise ValueError(f"Mixed polarities in manifest: {manifest_polarities}. Provide --polarity to override.")

        cfg = MassSightConfig(
            ppm=float(args.ppm),
            rt=float(args.rt),
            tight_ppm=float(args.tight_ppm),
            tight_rt=float(args.tight_rt),
            use_intensity=False,
            w_int=0.0,
            use_structure=False,
            polarity=polarity,
            mz_global_shift_model=str(getattr(args, "mz_global_shift_model", "auto")),
            mz_global_shift_anchor_ppm=float(getattr(args, "mz_global_shift_anchor_ppm", 5.0)),
            mz_global_shift_min_anchors=int(getattr(args, "mz_global_shift_min_anchors", 25)),
            mz_global_shift_min_gain=int(getattr(args, "mz_global_shift_min_gain", 5)),
            mz_global_shift_max_eval=int(getattr(args, "mz_global_shift_max_eval", 5000)),
            drift_fit_bootstrap=bool(args.drift_fit_bootstrap),
            drift_bootstrap_seed=int(args.drift_bootstrap_seed),
            rt_transfer_min_inlier_frac=float(getattr(args, "rt_transfer_min_inlier_frac", 0.02)),
            retention_nontransferable_policy=str(getattr(args, "retention_nontransferable_policy", "soft")),
            retention_nontransferable_weight=float(getattr(args, "retention_nontransferable_weight", 0.25)),
        )

        outputs, strategies = _run_cluster_strategies(
            studies=studies,
            cfg=cfg,
            out_dir=out_dir,
            strategy=str(args.strategy),
            hub_id=str(args.hub_id) if args.hub_id else None,
            hub_alignment=str(getattr(args, "hub_alignment", "auto")),
            full_join=bool(args.full_join),
            export_expr=bool(args.export_expr),
        )

        if bool(args.capture_manifest):
            ended_at = _utc_now_iso()
            runtime_seconds = time.perf_counter() - t0
            manifest_out = _write_run_manifest(
                args=args,
                started_at=started_at,
                ended_at=ended_at,
                runtime_seconds=runtime_seconds,
                outputs=outputs,
                extras={
                    "input_manifest": str(manifest_path),
                    "n_studies": int(len(studies)),
                    "strategies": [str(x) for x in strategies],
                },
            )
            outputs.append(str(manifest_out))

        print(f"Wrote clustering outputs to {out_dir}")
        return 0

    if args.cmd == "reuse":
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cache_dir = Path(args.cache_dir) if args.cache_dir else (out_dir / "cache")
        cache_dir.mkdir(parents=True, exist_ok=True)

        analysis_ids = _collect_analysis_ids(args)
        if not analysis_ids:
            raise ValueError("No analysis IDs provided. Use --analysis-id and/or --analysis-file.")

        bundles = []
        fetch_failures: list[dict[str, str]] = []
        for aid in analysis_ids:
            bundle = fetch_workbench_analysis(
                aid,
                cache_root=cache_dir,
                timeout_s=float(args.fetch_timeout_s),
                retries=int(args.fetch_retries),
                delay_s=float(args.fetch_delay_s),
                refresh=bool(args.refresh_fetch),
            )
            if bundle is None:
                fetch_failures.append({"analysis_id": str(aid), "reason": "fetch_failed_or_invalid_untarg"})
                continue
            bundles.append(bundle)

        if not bundles:
            raise RuntimeError("No analyses fetched successfully.")

        meta_rows = [dict(b.metadata) for b in bundles]
        metadata_df = pd.DataFrame(meta_rows).sort_values(["analysis_id"]).reset_index(drop=True)
        metadata_tsv = out_dir / "analysis_metadata.tsv"
        metadata_df.to_csv(metadata_tsv, sep="\t", index=False)

        grouped, dropped = group_workbench_analyses(
            bundles,
            min_group_size=int(args.min_group_size),
            allow_unknown_strata=bool(args.allow_unknown_strata),
        )
        dropped_all = dropped + fetch_failures

        outputs: list[str] = [str(metadata_tsv)]
        if dropped_all:
            dropped_df = pd.DataFrame(dropped_all).sort_values(["analysis_id", "reason"]).reset_index(drop=True)
            dropped_tsv = out_dir / "dropped_analyses.tsv"
            dropped_df.to_csv(dropped_tsv, sep="\t", index=False)
            outputs.append(str(dropped_tsv))

        groups_dir = out_dir / "groups"
        groups_dir.mkdir(parents=True, exist_ok=True)

        # Optional: export labeled metabolites (mwTab MS_METABOLITE_DATA, plus optional study_id targeted data).
        labeled_rows: list[dict[str, Any]] = []
        targeted_meta_by_study: dict[str, pd.DataFrame] = {}
        if bool(getattr(args, "export_labeled_metabolites", True)):
            # 1) mwTab MS_METABOLITE_DATA metabolites (always available when mwTab has MS_METABOLITE_DATA/Metabolites).
            for b in bundles:
                try:
                    mwtab_obj = json.loads(Path(str(b.mwtab_path)).read_text(encoding="utf-8"))
                except Exception:
                    mwtab_obj = None
                if not isinstance(mwtab_obj, dict):
                    continue
                rows = extract_mwtab_labeled_metabolites(mwtab_obj)
                for r in rows:
                    r2 = dict(r)
                    r2["analysis_id"] = str(b.analysis_id)
                    r2["study_id"] = str(b.metadata.get("study_id") or "")
                    labeled_rows.append(r2)

            # 2) Workbench study_id/.../data targeted matrices (optional fetch).
            if bool(getattr(args, "fetch_targeted_data", False)):
                max_bytes = int(float(getattr(args, "targeted_max_mib", 50.0)) * 1024 * 1024)
                study_ids = sorted(
                    {
                        str(b.metadata.get("study_id") or "").strip()
                        for b in bundles
                        if str(b.metadata.get("study_id") or "").strip()
                    }
                )
                for sid in study_ids:
                    path = fetch_workbench_study_data(
                        sid,
                        cache_root=cache_dir,
                        timeout_s=float(args.fetch_timeout_s),
                        retries=int(args.fetch_retries),
                        delay_s=float(args.fetch_delay_s),
                        refresh=bool(args.refresh_fetch),
                        max_bytes=max_bytes,
                    )
                    if path is None:
                        continue
                    try:
                        obj = json.loads(Path(str(path)).read_text(encoding="utf-8"))
                    except Exception:
                        continue
                    if not isinstance(obj, dict):
                        continue
                    meta_df, expr_df = parse_workbench_study_data_json(obj)
                    if meta_df.empty:
                        continue
                    targeted_meta_by_study[str(sid)] = meta_df

                    # Export per-study targeted matrices for downstream meta-analysis.
                    targeted_dir = out_dir / "targeted_data" / str(sid)
                    targeted_dir.mkdir(parents=True, exist_ok=True)
                    meta_path = targeted_dir / "named_metabolites.tsv"
                    meta_df.to_csv(meta_path, sep="\t", index=False)
                    outputs.append(str(meta_path))
                    if not expr_df.empty:
                        expr_path = targeted_dir / "named_metabolite_expr.tsv.gz"
                        expr_df.to_csv(expr_path, sep="\t", compression="gzip")
                        outputs.append(str(expr_path))

                    # Also record these as "labeled rows" keyed by metabolite_norm for overlap reporting.
                    for rec in meta_df.to_dict(orient="records"):
                        labeled_rows.append(dict(rec))

            # Write combined labeled metabolite table.
            labeled_tsv = out_dir / "analysis_labeled_metabolites.tsv"
            if labeled_rows:
                df = pd.DataFrame(labeled_rows)
                # Ensure consistent columns for consumers; leave unknown columns empty.
                want_cols = [
                    "analysis_id",
                    "study_id",
                    "source",
                    "metabolite_norm",
                    "metabolite_name",
                    "refmet_name",
                    "metabolite_id",
                    "pubchem_id",
                    "kegg_id",
                    "hmdb_id",
                    "chebi_id",
                    "inchikey",
                    "lipidmaps_id",
                    "other_id",
                    "other_id_type",
                    "units",
                    "analysis_summary",
                    "matrix_col",
                ]
                cols = [c for c in want_cols if c in df.columns] + [c for c in df.columns if c not in want_cols]
                df = df.loc[:, cols].copy()
                if "analysis_id" in df.columns:
                    df["analysis_id"] = df["analysis_id"].astype(str).map(lambda x: x.upper().strip())
                if "metabolite_norm" in df.columns:
                    df["metabolite_norm"] = df["metabolite_norm"].astype(str).map(normalize_metabolite_name)
                df = df.sort_values([c for c in ["analysis_id", "source", "metabolite_norm"] if c in df.columns]).reset_index(drop=True)
                df.to_csv(labeled_tsv, sep="\t", index=False)
                outputs.append(str(labeled_tsv))

        group_rows: list[dict[str, Any]] = []
        total_clustered = 0
        for group_key, members in sorted(grouped.items()):
            ion_mode_norm, _, chrom_norm = group_key.partition("__")
            polarity = str(args.polarity or ion_mode_norm).lower().strip()
            if polarity not in {"positive", "negative"}:
                # Keep behavior explicit; unknown polarity is skipped unless user overrides.
                for m in members:
                    dropped_all.append(
                        {
                            "analysis_id": str(m.analysis_id),
                            "group_key": str(group_key),
                            "reason": "unknown_polarity_requires_override",
                        }
                    )
                continue

            group_out = groups_dir / group_key
            group_out.mkdir(parents=True, exist_ok=True)
            group_manifest_path = group_out / "group_manifest.json"

            specs = [bundle_to_study_spec(m, polarity_override=polarity) for m in members]
            group_manifest_obj = {
                "schema_version": "mass_sight.reuse_group_manifest.v1",
                "group_key": str(group_key),
                "ion_mode_norm": str(ion_mode_norm),
                "chromatography_norm": str(chrom_norm),
                "studies": [
                    {
                        "analysis_id": str(s.analysis_id),
                        "feature_format": str(s.feature_format),
                        "matrix_path": str(s.matrix_path),
                        "polarity": str(s.polarity),
                        "chromatography": str(s.chromatography or ""),
                        "metadata": dict(s.metadata),
                    }
                    for s in specs
                ],
            }
            group_manifest_path.write_text(json.dumps(group_manifest_obj, indent=2, sort_keys=True) + "\n")
            outputs.append(str(group_manifest_path))

            # Shared labeled metabolites (within-group) for downstream targeted/meta-analysis.
            if bool(getattr(args, "export_labeled_metabolites", True)) and labeled_rows:
                analysis_ids_in_group = [str(m.analysis_id) for m in members]
                # Build per-analysis best row per metabolite_norm with a simple source priority.
                pri = {"workbench_study_data": 0, "mwtab_ms_metabolite_data": 1}
                rows_by_an_norm: dict[tuple[str, str], dict[str, Any]] = {}
                for rec in labeled_rows:
                    an = str(rec.get("analysis_id") or "").upper().strip()
                    if an not in analysis_ids_in_group:
                        continue
                    norms = {
                        normalize_metabolite_name(rec.get("metabolite_norm") or ""),
                        normalize_metabolite_name(rec.get("refmet_name") or ""),
                        normalize_metabolite_name(rec.get("metabolite_name") or ""),
                    }
                    norms = {n for n in norms if n}
                    if not norms:
                        continue
                    for norm in sorted(norms):
                        key = (an, norm)
                        prev = rows_by_an_norm.get(key)
                        if prev is None:
                            rows_by_an_norm[key] = dict(rec)
                            continue
                        p_prev = pri.get(str(prev.get("source") or ""), 99)
                        p_new = pri.get(str(rec.get("source") or ""), 99)
                        if p_new < p_prev:
                            rows_by_an_norm[key] = dict(rec)

                norms_by_an: dict[str, set[str]] = {}
                for an in analysis_ids_in_group:
                    norms_by_an[an] = {norm for (a, norm) in rows_by_an_norm.keys() if a == an}

                shared_dir = group_out / "shared_labeled_metabolites"
                shared_dir.mkdir(parents=True, exist_ok=True)

                # Intersection across all studies in the group.
                all_shared = set.intersection(*[norms_by_an[an] for an in analysis_ids_in_group]) if analysis_ids_in_group else set()
                if all_shared:
                    out_all = []
                    for norm in sorted(all_shared):
                        # Use the first study's row as representative.
                        rec0 = rows_by_an_norm.get((analysis_ids_in_group[0], norm), {})
                        out_all.append(
                            {
                                "metabolite_norm": norm,
                                "metabolite_name": str(rec0.get("refmet_name") or rec0.get("metabolite_name") or ""),
                            }
                        )
                    all_path = shared_dir / "shared_all.tsv"
                    pd.DataFrame(out_all).to_csv(all_path, sep="\t", index=False)
                    outputs.append(str(all_path))

                # Pairwise shared metabolite lists (one file per pair) + summary counts.
                pair_rows = []
                for i, an1 in enumerate(sorted(analysis_ids_in_group)):
                    for an2 in sorted(analysis_ids_in_group)[i + 1 :]:
                        shared = sorted(norms_by_an.get(an1, set()) & norms_by_an.get(an2, set()))
                        pair_rows.append({"analysis_id_a": an1, "analysis_id_b": an2, "n_shared": int(len(shared))})
                        if not shared:
                            continue
                        out = []
                        for norm in shared:
                            r1 = rows_by_an_norm.get((an1, norm), {})
                            r2 = rows_by_an_norm.get((an2, norm), {})
                            out.append(
                                {
                                    "metabolite_norm": norm,
                                    "name_a": str(r1.get("refmet_name") or r1.get("metabolite_name") or ""),
                                    "name_b": str(r2.get("refmet_name") or r2.get("metabolite_name") or ""),
                                    "source_a": str(r1.get("source") or ""),
                                    "source_b": str(r2.get("source") or ""),
                                    "metabolite_id_a": str(r1.get("metabolite_id") or ""),
                                    "metabolite_id_b": str(r2.get("metabolite_id") or ""),
                                    "pubchem_id_a": str(r1.get("pubchem_id") or ""),
                                    "pubchem_id_b": str(r2.get("pubchem_id") or ""),
                                    "kegg_id_a": str(r1.get("kegg_id") or ""),
                                    "kegg_id_b": str(r2.get("kegg_id") or ""),
                                    "inchikey_a": str(r1.get("inchikey") or ""),
                                    "inchikey_b": str(r2.get("inchikey") or ""),
                                }
                            )
                        pair_path = shared_dir / f"{an1}__{an2}.tsv"
                        pd.DataFrame(out).to_csv(pair_path, sep="\t", index=False)
                        outputs.append(str(pair_path))

                if pair_rows:
                    pairs_path = shared_dir / "shared_pair_counts.tsv"
                    pd.DataFrame(pair_rows).sort_values(["analysis_id_a", "analysis_id_b"]).to_csv(pairs_path, sep="\t", index=False)
                    outputs.append(str(pairs_path))

            try:
                studies = [load_study(s) for s in specs]
                use_intensity, intensity_standardize, intensity_reason = resolve_use_intensity_mode(
                    [m.metadata for m in members],
                    mode=str(args.use_intensity_mode),
                )

                cfg = MassSightConfig(
                    ppm=float(args.ppm),
                    rt=float(args.rt),
                    tight_ppm=float(args.tight_ppm),
                    tight_rt=float(args.tight_rt),
                    use_intensity=bool(use_intensity),
                    w_int=1.0 if bool(use_intensity) else 0.0,
                    intensity_standardize=str(intensity_standardize),
                    use_structure=False,
                    polarity=str(polarity),
                    mz_global_shift_model=str(getattr(args, "mz_global_shift_model", "auto")),
                    mz_global_shift_anchor_ppm=float(getattr(args, "mz_global_shift_anchor_ppm", 5.0)),
                    mz_global_shift_min_anchors=int(getattr(args, "mz_global_shift_min_anchors", 25)),
                    mz_global_shift_min_gain=int(getattr(args, "mz_global_shift_min_gain", 5)),
                    mz_global_shift_max_eval=int(getattr(args, "mz_global_shift_max_eval", 5000)),
                    drift_fit_bootstrap=bool(args.drift_fit_bootstrap),
                    drift_bootstrap_seed=int(args.drift_bootstrap_seed),
                )

                group_cluster_outputs, strategies = _run_cluster_strategies(
                    studies=studies,
                    cfg=cfg,
                    out_dir=group_out,
                    strategy=str(args.strategy),
                    hub_id=None,
                    hub_alignment=str(args.hub_alignment),
                    full_join=bool(args.full_join),
                    export_expr=bool(args.export_expr),
                )
                outputs.extend(group_cluster_outputs)
                total_clustered += int(len(studies))

                group_rows.append(
                    {
                        "group_key": str(group_key),
                        "n_studies": int(len(studies)),
                        "ion_mode_norm": str(ion_mode_norm),
                        "chromatography_norm": str(chrom_norm),
                        "use_intensity_requested": str(args.use_intensity_mode),
                        "use_intensity_effective": str("enabled" if use_intensity else "disabled"),
                        "use_intensity_reason": str(intensity_reason),
                        "strategies": "|".join(strategies),
                        "group_out_dir": str(group_out),
                    }
                )
            except Exception as exc:
                for m in members:
                    dropped_all.append(
                        {
                            "analysis_id": str(m.analysis_id),
                            "group_key": str(group_key),
                            "reason": "group_processing_failed",
                            "error": str(exc),
                        }
                    )

        if dropped_all:
            dropped_df = pd.DataFrame(dropped_all).sort_values(["analysis_id", "reason"]).reset_index(drop=True)
            dropped_tsv = out_dir / "dropped_analyses.tsv"
            dropped_df.to_csv(dropped_tsv, sep="\t", index=False)
            if str(dropped_tsv) not in outputs:
                outputs.append(str(dropped_tsv))

        groups_tsv = out_dir / "reuse_groups.tsv"
        groups_df = pd.DataFrame(
            group_rows,
            columns=[
                "group_key",
                "n_studies",
                "ion_mode_norm",
                "chromatography_norm",
                "use_intensity_requested",
                "use_intensity_effective",
                "use_intensity_reason",
                "strategies",
                "group_out_dir",
            ],
        )
        groups_df.to_csv(groups_tsv, sep="\t", index=False)
        outputs.append(str(groups_tsv))

        summary = {
            "schema_version": "mass_sight.reuse_summary.v1",
            "n_requested_analyses": int(len(analysis_ids)),
            "n_fetched_analyses": int(len(bundles)),
            "n_clustered_analyses": int(total_clustered),
            "n_groups": int(len(group_rows)),
            "n_dropped": int(len(dropped_all)),
            "groups": group_rows,
        }
        summary_path = out_dir / "reuse_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        outputs.append(str(summary_path))

        if bool(args.capture_manifest):
            ended_at = _utc_now_iso()
            runtime_seconds = time.perf_counter() - t0
            manifest_out = _write_run_manifest(
                args=args,
                started_at=started_at,
                ended_at=ended_at,
                runtime_seconds=runtime_seconds,
                outputs=outputs,
                extras={
                    "n_requested_analyses": int(len(analysis_ids)),
                    "n_fetched_analyses": int(len(bundles)),
                    "n_groups": int(len(group_rows)),
                },
            )
            outputs.append(str(manifest_out))

        print(f"Wrote reuse outputs to {out_dir}")
        return 0

    if args.cmd == "find":
        out_path = Path(str(args.out))
        out_path.parent.mkdir(parents=True, exist_ok=True)
        cache_dir = Path(args.cache_dir) if args.cache_dir else (out_path.parent / "cache")
        cache_dir.mkdir(parents=True, exist_ok=True)

        # Lazy import so that the rest of the CLI stays usable even if TUI deps are missing.
        try:
            from .workbench_finder_tui import run_workbench_finder_tui
        except Exception as exc:  # pragma: no cover
            raise RuntimeError(
                "Textual TUI dependencies are missing or failed to import. "
                "Install the `textual` package to use `mass_sight find`."
            ) from exc

        saved_path = run_workbench_finder_tui(
            out_path=out_path,
            cache_dir=cache_dir,
            refresh=bool(args.refresh_fetch),
            timeout_s=float(args.fetch_timeout_s),
            retries=int(args.fetch_retries),
            delay_s=float(args.fetch_delay_s),
        )
        if saved_path is not None:
            print(f"Wrote Workbench selection manifest to {saved_path}")
        else:
            print("Exited without saving selection manifest (press 's' from an analyses screen to save).")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
