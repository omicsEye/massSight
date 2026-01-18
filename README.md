# massSight

massSight implements probabilistic drift-aware optimal transport for cross-study LC–MS feature matching.

[![](https://zenodo.org/badge/608216683.svg)](https://doi.org/10.5281/zenodo.8101763)

## Install

Requires Python ≥3.10.

After installing, you can:

- import `mass_sight` in Python
- run `mass_sight` on the command line

### Using `uv`

```bash
uv pip install git+https://github.com/omicsEye/massSight
```

### Using `pip`

```bash
pip install git+https://github.com/omicsEye/massSight
```

Core dependencies include `numpy`, `pandas`, `scikit-learn`, and `scipy`.

## Quickstart

### Input feature tables

`massSight` expects two feature tables with at least:

- `MZ` (m/z)
- `RT` (retention time in minutes)
- `Intensity` (a per-feature intensity summary; optional)

### `MassSightConfig` and column names

`massSight` is configured via `MassSightConfig`. By default it expects canonical column names (`MZ`, `RT`, `Intensity`).

If your tables use different column names, either rename your columns or pass a schema override:

```python
from mass_sight import MassSightConfig

cfg = MassSightConfig(mz_col="mz", rt_col="rt_min", intensity_col="area")
```

If `study_a` and `study_b` use different schemas, use study-specific overrides:

```python
cfg = MassSightConfig(
    mz_col_study_a="mz",
    rt_col_study_a="rt_min",
    mz_col_study_b="m_over_z",
    rt_col_study_b="rt",
)
```

For the CLI, use `--mz-col`, `--rt-col`, and `--intensity-col`.

### Python usage

```python
import pandas as pd
from mass_sight import MassSightConfig, match_features

study_a = pd.read_csv("study_a.csv")
study_b = pd.read_csv("study_b.csv")

cfg = MassSightConfig()
res = match_features(study_a, study_b, cfg)

top1 = res.top1 # id1, id2, decision, support, margin, prob_raw, ...
candidates = res.candidates  # residuals, log-likelihoods, OT weights, etc.
```

### Command-line usage

```bash
mass_sight match study_a.csv study_b.csv \
  --out-candidates candidates.csv \
  --out-top1 top1.csv
```

## Citation

- See `CITATION.cff`.

## License

MIT.
