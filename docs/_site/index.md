# massSight


This package provides:

- Supervised matching with boosted trees (LightGBM), calibrated with
  out‑of‑fold isotonic regression \[1\]
- Per‑row probability mass estimates for matches and non-matches
  (suitable for sampling/MI).
- MI fold‑change with randomized global 1–1 matchings using
  Gumbel‑perturbed Hungarian \[2\]

## Install

``` bash
pip install git+https://github.com/omicsEye/massSight.git
```

Requires Python $\geq$ 3.10 and the following core deps: `numpy`,
`pandas`, `scikit-learn`, `scipy`, `lightgbm`.

## Quickstart (Python API)

### 1. Matching two feature tables

``` python
import pandas as pd
from mass_sight import match_ml, MLMatchConfig

ds1 = pd.read_csv("DS1.csv") 
ds2 = pd.read_csv("DS2.csv")

cfg = MLMatchConfig(
    ppm=12.0, 
    rt=1.0, 
    model="auto", 
    seed=123,
    mz_col="MZ", 
    rt_col="RT", 
    annotation_col="Annotation_ID", 
    compound_id_col="Compound_ID",
    # Provide either one intensity col, 
    # or a list/regex to average (then log10):
    intensity_col=None,
    intensity_cols=["S1_int", "S2_int", "S3_int"],  # example
    intensity_regex=None,
)

res = match_ml(ds1, ds2, cfg)
cand = res.candidates  # DataFrame: id1, id2, prob_cal, p_row, p0, ppm_diff, rt_diff, roi_diff, log_int_diff

# Optional 1–1 predictions
top1 = res.top1         # greedy best per row
hung = res.hungarian    # global 1–1 (Hungarian)
```

Intensity handling:

- If `intensity_col` is set, Intensity_log10 = log10(intensity + eps).
- If `intensity_cols` or `intensity_regex` is set, Intensity_log10 =
  log10(mean(raw intensities) + eps).
- Features use ppm_diff, rt_diff, retention‑order diff (Δro), and
  log_intensity diff (Δ log10 mean).

### 2. MI fold‑change with per‑sample matrices

``` python
from mass_sight import mi_fold_change

ds1_expr = pd.read_csv("DS1_expr.csv", index_col=0)  # rows=DS row indices (or map beforehand), cols=sample IDs
ds1_meta = pd.read_csv("DS1_samples.csv")           # columns: sample, group (2 levels)
ds2_expr = pd.read_csv("DS2_expr.csv", index_col=0)
ds2_meta = pd.read_csv("DS2_samples.csv")

mi = mi_fold_change(ds1_expr, ds1_meta, ds2_expr, ds2_meta, cand, M=200, seed=42)
# mi columns: entity_id, id1, id2, logFC_hat, SE_MI, CI_low, CI_high, FC, draws_used, match_prob, p_row_best, margin
```

What MI does:

- Samples per‑row match vs no‑match from p_row/p0.
- For sampled match rows, solves a randomized global 1–1 matching via
  Hungarian on −log(p) + Gumbel noise.
- Fits OLS on concatenated samples: log1p(intensity) ~ 1 + group +
  dataset.
- Combines M draws via Rubin’s rules to produce logFC_hat, SE, 95% CI,
  and FC.

## CLI (optional)

After install, you can run:

``` bash
# 1. Matching with schema mapping and intensity rules
mass-sight match DS1.csv DS2.csv \
  --out candidates.csv \
  --ppm 12 --rt 1.0 \
  --mz-col MZ --rt-col RT \ 
  --annotation-col Annotation_ID \
  --compound-id-col Compound_ID \
  --int-col Intensity # or: --int-col S1_int --int-col S2_int  ...
  # or: --int-regex '^S\\d+_int$'

# 2. MI fold-change
mass-sight mi \
  --ds1-expr DS1_expr.csv --ds1-meta DS1_samples.csv \
  --ds2-expr DS2_expr.csv --ds2-meta DS2_samples.csv \
  --candidates candidates.csv --M 200 --out mi_fc.csv
```

## Notes & best practices

- Matching windows (ppm, rt) are the only knobs you usually need.
- Provide intensity columns or regex to compute a stable log10 mean per
  feature; the model uses Δlog10 mean as one of the features.
- For MI, ensure both datasets have group labels with two levels and
  adequate replication.
- The dataset offset in the MI model absorbs platform/global shifts; do
  not pre‑normalize across datasets.

## License & citation

MIT. If you use massSight, please cite the repo and version tag you
used.

## References

<div id="refs" class="references csl-bib-body" entry-spacing="0">

<div id="ref-zadroznyTransformingClassifierScores2002"
class="csl-entry">

<span class="csl-left-margin">\[1\]
</span><span class="csl-right-inline">B. Zadrozny and C. Elkan,
“Transforming classifier scores into accurate multiclass probability
estimates,” in *Proceedings of the eighth ACM SIGKDD international
conference on Knowledge discovery and data mining*, in KDD ’02. New
York, NY, USA: Association for Computing Machinery, Jul. 2002, pp.
694–699. doi:
[10.1145/775047.775151](https://doi.org/10.1145/775047.775151).</span>

</div>

<div id="ref-liEfficientFeatureLearning" class="csl-entry">

<span class="csl-left-margin">\[2\]
</span><span class="csl-right-inline">K. Li, K. Swersky, and R. Zemel,
“Efficient <span class="nocase">Feature Learning Using
Perturb-and-MAP</span>.”</span>

</div>

</div>
