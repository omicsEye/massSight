import argparse
import sys
from pathlib import Path
import pandas as pd
from .ml_matcher import match_ml, MLMatchConfig
from .multi_match import align_multi, MultiMatchConfig
from .mi_fold_change import mi_fold_change


def _add_match_parser(sub):
    p = sub.add_parser("match", help="Match two datasets with ML + calibration")
    p.add_argument("ds1", type=str, help="Path to DS1 CSV")
    p.add_argument("ds2", type=str, help="Path to DS2 CSV")
    p.add_argument("--out", type=str, required=True, help="Output CSV for calibrated candidates")
    p.add_argument("--ppm", type=float, default=12.0)
    p.add_argument("--rt", type=float, default=1.0)
    # schema
    p.add_argument("--mz-col", type=str, default="MZ")
    p.add_argument("--rt-col", type=str, default="RT")
    p.add_argument("--annotation-col", type=str, default="Annotation_ID")
    p.add_argument("--compound-id-col", type=str, default="Compound_ID")
    # intensity options
    p.add_argument("--int-col", type=str, default=None, help="Single intensity column name")
    p.add_argument("--int-cols", dest="int_cols", action="append", default=None, help="Repeat for multiple columns (mean log10)")
    p.add_argument("--int-regex", type=str, default=None, help="Regex to select multiple intensity columns")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--save-model", dest="save_model", type=str, default=None, help="Path to save fitted model artifact (joblib)")
    return p


def _add_mi_parser(sub):
    p = sub.add_parser("mi", help="Multiple imputation fold-change on matched candidates")
    p.add_argument("--ds1-expr", required=True, type=str, help="DS1 expression matrix CSV (rows=features, cols=samples)")
    p.add_argument("--ds1-meta", required=True, type=str, help="DS1 sample metadata CSV (columns: sample, group)")
    p.add_argument("--ds2-expr", required=True, type=str, help="DS2 expression matrix CSV")
    p.add_argument("--ds2-meta", required=True, type=str, help="DS2 sample metadata CSV")
    p.add_argument("--candidates", required=True, type=str, help="Calibrated candidates CSV from match")
    p.add_argument("--out", required=True, type=str, help="Output CSV for MI fold-change results")
    p.add_argument("--M", type=int, default=100)
    p.add_argument("--seed", type=int, default=123)
    return p


def main(argv=None):
    argv = argv or sys.argv[1:]
    ap = argparse.ArgumentParser(prog="mass-sight", description="ML matching + MI fold-change for LC–MS")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_match = _add_match_parser(sub)
    p_mi = _add_mi_parser(sub)
    p_multi = _add_multi_parser(sub)
    args = ap.parse_args(argv)

    if args.cmd == "match":
        ds1 = pd.read_csv(args.ds1)
        ds2 = pd.read_csv(args.ds2)
        cfg = MLMatchConfig(
            ppm=args.ppm, rt=args.rt, seed=args.seed,
            mz_col=args.mz_col, rt_col=args.rt_col,
            annotation_col=args.annotation_col, compound_id_col=args.compound_id_col,
            intensity_col=None if args.int_cols else args.int_col,
            intensity_cols=args.int_cols, intensity_regex=args.int_regex,
            save_model_path=args.save_model,
        )
        res = match_ml(ds1, ds2, cfg)
        out = res.candidates
        out.to_csv(args.out, index=False)
        print(f"wrote {args.out} with {len(out)} candidate edges")
        return 0

    if args.cmd == "mi":
        ds1_expr = pd.read_csv(args.ds1_expr, index_col=0)
        ds1_meta = pd.read_csv(args.ds1_meta)
        ds2_expr = pd.read_csv(args.ds2_expr, index_col=0)
        ds2_meta = pd.read_csv(args.ds2_meta)
        cand = pd.read_csv(args.candidates)
        res = mi_fold_change(ds1_expr, ds1_meta, ds2_expr, ds2_meta, cand, M=args.M, seed=args.seed)
        res.to_csv(args.out, index=False)
        print(f"wrote {args.out} with {len(res)} entities")
        return 0

    if args.cmd == "multi":
        if len(args.datasets) < 3:
            print("Need at least 3 datasets for multi alignment", file=sys.stderr)
            return 2
        dfs = [pd.read_csv(p) for p in args.datasets]
        names = [Path(p).stem if "/" in p or "." in p else str(p) for p in args.datasets]
        mcfg = MultiMatchConfig(
            ppm=args.ppm, rt=args.rt, seed=args.seed,
            topk=args.topk, tau_main=args.tau, prob_floor=args.prob_floor, max_tuples_per_anchor=args.max_tuples,
            mz_col=args.mz_col, rt_col=args.rt_col, annotation_col=args.annotation_col, compound_id_col=args.compound_id_col,
            intensity_col=None if args.int_cols else args.int_col, intensity_cols=args.int_cols, intensity_regex=args.int_regex,
        )
        res = align_multi(dfs, names=names, config=mcfg)
        out = res.clusters
        out.to_csv(args.out, index=False)
        print(f"wrote {args.out} with {len(out)} clusters; K={len(dfs)} topk={args.topk}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
def _add_multi_parser(sub):
    p = sub.add_parser("multi", help="Align K>=3 datasets globally using calibrated pairwise ML")
    p.add_argument("datasets", nargs="+", help="CSV paths for datasets in order (e.g., DS1 DS2 DS3 DS4)")
    p.add_argument("--out", required=True, help="Output CSV for selected K-way clusters")
    p.add_argument("--ppm", type=float, default=12.0)
    p.add_argument("--rt", type=float, default=1.0)
    p.add_argument("--topk", type=int, default=3)
    p.add_argument("--tau", type=float, default=0.2)
    p.add_argument("--prob-floor", dest="prob_floor", type=float, default=0.05)
    p.add_argument("--max-tuples", dest="max_tuples", type=int, default=50)
    # schema (shared)
    p.add_argument("--mz-col", type=str, default="MZ")
    p.add_argument("--rt-col", type=str, default="RT")
    p.add_argument("--annotation-col", type=str, default="Annotation_ID")
    p.add_argument("--compound-id-col", type=str, default="Compound_ID")
    p.add_argument("--int-col", type=str, default=None)
    p.add_argument("--int-cols", dest="int_cols", action="append", default=None)
    p.add_argument("--int-regex", type=str, default=None)
    p.add_argument("--seed", type=int, default=123)
    return p
