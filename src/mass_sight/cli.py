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
    group_workbench_analyses,
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
        help="Optional text file with one analysis ID per line.",
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


def _collect_analysis_ids(args: argparse.Namespace) -> list[str]:
    values: list[str] = []
    values.extend([str(x).strip() for x in (args.analysis_ids or []) if str(x).strip()])
    if getattr(args, "analysis_file", None):
        path = Path(str(args.analysis_file))
        if not path.exists():
            raise FileNotFoundError(path)
        for line in path.read_text(encoding="utf-8").splitlines():
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
            drift_fit_bootstrap=bool(args.drift_fit_bootstrap),
            drift_bootstrap_seed=int(args.drift_bootstrap_seed),
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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
