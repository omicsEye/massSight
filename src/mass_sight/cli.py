import argparse
import sys
import pandas as pd
from .ml_matcher import match_mlnodrift, MLMatchConfig
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
    p.add_argument("--int-col", dest="int_cols", action="append", default=None, help="Repeat for multiple columns (mean log10)")
    p.add_argument("--int-regex", type=str, default=None, help="Regex to select multiple intensity columns")
    p.add_argument("--seed", type=int, default=123)
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
        )
        res = match_mlnodrift(ds1, ds2, cfg)
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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())

