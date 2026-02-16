# Reviewer demo (runnable)

This bundle is a minimal, deterministic smoke test for the `mass_sight` CLI.

## Run

```bash
uv run mass_sight match \
  examples/reviewer_demo/study_a.csv \
  examples/reviewer_demo/study_b.csv \
  --out-candidates examples/reviewer_demo/output/candidates.csv \
  --out-top1 examples/reviewer_demo/output/top1.csv
```

Default CLI behavior also writes a run manifest:

- `examples/reviewer_demo/output/candidates.run_manifest.json`

## Verify expected outputs

```bash
uv run python examples/reviewer_demo/verify_outputs.py \
  --candidates examples/reviewer_demo/output/candidates.csv \
  --top1 examples/reviewer_demo/output/top1.csv
```

Expected summary counts are defined in:

- `examples/reviewer_demo/expected_summary.json`
