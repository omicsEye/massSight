#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    p = argparse.ArgumentParser(description="Verify reviewer demo outputs against expected summary counts.")
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--top1", type=Path, required=True)
    p.add_argument("--expected", type=Path, default=Path(__file__).with_name("expected_summary.json"))
    args = p.parse_args()

    expected = json.loads(args.expected.read_text())
    cand = pd.read_csv(args.candidates)
    top1 = pd.read_csv(args.top1)

    got = {
        "n_candidates": int(len(cand)),
        "n_top1_rows": int(len(top1)),
        "n_top1_matches": int((top1["decision"].astype(str) == "match").sum()) if "decision" in top1.columns else 0,
        "n_top1_no_candidate": int((top1["decision"].astype(str) == "no_candidate").sum()) if "decision" in top1.columns else 0,
    }
    for key, value in got.items():
        if int(expected.get(key, -1)) != int(value):
            raise SystemExit(f"Mismatch for {key}: expected={expected.get(key)} got={value}")

    print("Reviewer demo outputs verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
