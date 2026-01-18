import argparse
import sys

import pandas as pd

from .matcher import MassSightConfig, match_features
from .multistudy import (
    cluster_hub_barycenter_ot,
    cluster_hub_mutual_top1,
    cluster_symmetric_mutual_graph,
    load_manifest,
    load_study,
 )


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

    # Candidate expansion (ion-aware mass offsets)
    p.add_argument("--mz-shift-mode", choices=["none", "auto", "manual"], default="auto")
    p.add_argument(
        "--mz-shift-delta-da",
        dest="mz_shift_deltas_da",
        action="append",
        default=None,
        help="Manual discrete m/z shift in Da (repeatable; used when --mz-shift-mode manual)",
    )
    p.add_argument("--mz-shift-penalty", type=float, default=0.5, help="Penalty for non-zero shift hypotheses")

    # Intensity / unmatched mass
    p.add_argument("--use-intensity", dest="use_intensity", action="store_true", help="Enable intensity residual term")
    p.add_argument("--no-intensity", dest="use_intensity", action="store_false",
                   help="Disable intensity residual term (default)")
    p.set_defaults(use_intensity=False)
    p.add_argument(
        "--consensus",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Output a consensus 1–1 alignment (default: disabled; use OT top‑1 + explicit no-match).",
    )
    p.add_argument(
        "--consensus-samples",
        type=int,
        default=50,
        help="Number of sampled 1–1 matchings for consensus stability (default: 50).",
    )
    p.add_argument(
        "--consensus-min-support",
        type=float,
        default=0.25,
        help="Consensus output: minimum edge inclusion frequency to accept a match (default: 0.25).",
    )
    p.add_argument(
        "--consensus-gumbel-scale",
        type=float,
        default=1.0,
        help="Gumbel noise scale (temperature) for consensus sampling (default: 1.0).",
    )
    p.add_argument(
        "--consensus-seed",
        type=int,
        default=0,
        help="RNG seed for consensus sampling (default: 0).",
    )
    p.add_argument(
        "--allow-unmatched",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add a null column to allow 'no match' (default: enabled)",
    )
    p.add_argument(
        "--abstain",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable discriminative abstention (support + OT margin) (default: disabled).",
    )
    p.add_argument(
        "--abstain-support",
        type=float,
        default=0.25,
        help="Abstain when support < this AND margin < --abstain-margin (default: 0.25).",
    )
    p.add_argument(
        "--abstain-margin",
        type=float,
        default=0.1,
        help="Abstain when support < --abstain-support AND margin < this (default: 0.1).",
    )
    p.add_argument("--null-mass", type=float, default=0.1, help="OT mass assigned to the null column when enabled")

    # Schema
    p.add_argument("--mz-col", type=str, default="MZ")
    p.add_argument("--rt-col", type=str, default="RT")
    p.add_argument("--intensity-col", type=str, default="Intensity")
    p.add_argument("--polarity", choices=["positive", "negative"], default="positive")
    return p


def _add_cluster_parser(sub):
    p = sub.add_parser("cluster", help="Cluster multiple untargeted studies with massSight")
    p.add_argument("--manifest", required=True, type=str, help="Path to YAML/JSON multi-study manifest")
    p.add_argument("--out-dir", required=True, type=str, help="Output directory")
    p.add_argument(
        "--cluster-consensus",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use consensus 1–1 output within pairwise matching (default: disabled).",
    )
    p.add_argument("--cluster-consensus-samples", type=int, default=50)
    p.add_argument("--cluster-consensus-min-support", type=float, default=0.25)
    p.add_argument("--cluster-consensus-gumbel-scale", type=float, default=1.0)
    p.add_argument("--cluster-consensus-seed", type=int, default=0)
    p.add_argument(
        "--cluster-abstain",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable discriminative abstention within pairwise matching (default: disabled).",
    )
    p.add_argument("--cluster-abstain-support", type=float, default=0.25)
    p.add_argument("--cluster-abstain-margin", type=float, default=0.1)

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
                   help="Override polarity for shift hypotheses (default: infer/validate from manifest).")

    p.add_argument("--drift-fit-bootstrap", action="store_true", default=False,
                   help="Fit drift models on a bootstrap resample of tight matches (perturb-and-match).")
    p.add_argument("--drift-bootstrap-seed", type=int, default=0, help="RNG seed for drift bootstrap resampling.")
    return p


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    ap = argparse.ArgumentParser(prog="mass_sight", description="massSight matching for LC–MS")
    sub = ap.add_subparsers(dest="cmd", required=True)
    _add_match_parser(sub)
    _add_cluster_parser(sub)
    args = ap.parse_args(argv)

    if args.cmd == "match":
        study_a = pd.read_csv(args.study_a)
        study_b = pd.read_csv(args.study_b)
        cfg = MassSightConfig(
            ppm=float(args.ppm),
            rt=float(args.rt),
            tight_ppm=float(args.tight_ppm),
            tight_rt=float(args.tight_rt),
            mz_shift_mode=str(args.mz_shift_mode),
            mz_shift_manual_deltas_da=[float(x) for x in (args.mz_shift_deltas_da or [])],
            mz_shift_penalty=float(args.mz_shift_penalty),
            use_intensity=bool(args.use_intensity),
            use_consensus=bool(args.consensus),
            consensus_samples=int(args.consensus_samples),
            consensus_min_support=float(args.consensus_min_support),
            consensus_gumbel_scale=float(args.consensus_gumbel_scale),
            consensus_seed=int(args.consensus_seed),
            use_discriminative_abstention=bool(args.abstain),
            abstain_support=float(args.abstain_support),
            abstain_margin=float(args.abstain_margin),
            allow_unmatched=bool(args.allow_unmatched),
            null_mass=float(args.null_mass),
            mz_col=str(args.mz_col),
            rt_col=str(args.rt_col),
            intensity_col=str(args.intensity_col) if args.intensity_col else None,
            polarity=str(args.polarity),
        )
        res = match_features(study_a, study_b, cfg)
        res.candidates.to_csv(args.out_candidates, index=False)
        if args.out_top1:
            res.top1.to_csv(args.out_top1, index=False)
        print(
            f"wrote {args.out_candidates} ({len(res.candidates)} candidates); "
            f"top1 rows={len(res.top1)}"
        )
        return 0

    if args.cmd == "cluster":
        from pathlib import Path
        import json

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
            use_consensus=bool(args.cluster_consensus),
            consensus_samples=int(args.cluster_consensus_samples),
            consensus_min_support=float(args.cluster_consensus_min_support),
            consensus_gumbel_scale=float(args.cluster_consensus_gumbel_scale),
            consensus_seed=int(args.cluster_consensus_seed),
            use_discriminative_abstention=bool(args.cluster_abstain),
            abstain_support=float(args.cluster_abstain_support),
            abstain_margin=float(args.cluster_abstain_margin),
        )

        strategies = [args.strategy] if args.strategy != "both" else ["hub", "symmetric"]
        for strat in strategies:
            if strat == "hub":
                hub_id = str(args.hub_id) if args.hub_id else None
                hub_alignment = str(getattr(args, "hub_alignment", "auto")).lower().strip()
                use_barycenter = hub_alignment == "barycenter" or (
                    hub_alignment == "auto" and len(studies) >= 3
                )

                if use_barycenter:
                    res = cluster_hub_barycenter_ot(
                        studies,
                        cfg,
                        hub_id=hub_id,
                        include_unmatched=bool(args.full_join),
                        export_expr=bool(args.export_expr),
                    )
                else:
                    res = cluster_hub_mutual_top1(
                        studies,
                        cfg,
                        hub_id=hub_id,
                        export_expr=bool(args.export_expr),
                        include_unmatched=bool(args.full_join),
                    )
            else:
                res = cluster_symmetric_mutual_graph(studies, cfg, export_expr=bool(args.export_expr))

            sub_out = out_dir / res.strategy if args.strategy == "both" else out_dir
            sub_out.mkdir(parents=True, exist_ok=True)

            res.cluster_metadata.to_csv(sub_out / "cluster_metadata.tsv", sep="\t", index=False)
            res.cluster_map.to_csv(sub_out / "cluster_map.tsv", sep="\t", index=False)
            for study_id, mat in res.cluster_expr.items():
                mat.to_csv(sub_out / f"cluster_expr_{study_id}.tsv.gz", sep="\t", compression="gzip")
            (sub_out / "diagnostics.json").write_text(json.dumps(res.diagnostics, indent=2, sort_keys=True))

        print(f"Wrote clustering outputs to {out_dir}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
