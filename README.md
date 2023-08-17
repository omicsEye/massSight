# massSight

- [Examples](#examples)
- [Description](#description)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
  - [The `massSight` Object](#ms-obj)
- [Alignment](#align)
  - [`auto_combine()`](#auto_combine)
  - [`ml_match()`](#ml_match)
  - [Plotting results from alignment](#plotting-results-from-alignment)
- [Dev Instructions](#dev-instrutions)
  - [Installation](#installation-1)

<!-- README.md is generated from README.qmd. Please edit that file -->

<img src="man/figures/massSight.png" align="right" width="30%"/></a>

<div>

[![](https://zenodo.org/badge/608216683.svg)](https://zenodo.org/badge/latestdoi/608216683)

DOI

</div>

`massSight` is an R package for combining and scaling of LC-MS
metabolomics data.

- Citation: if you use `massSight`, please cite our manuscript: Chiraag
  Gohel and Ali Rahnavard. (2023). massSight: Metabolomics meta-analysis
  through multi-study data scaling, integration, and harmonization.
  <https://github.com/omicsEye/massSight>

## Examples

Examples and extensive documentation can be found
[here](omicseye.github.io/massSight/)

## Description

## Installation

``` r
devtools::install_github("omicsEye/massSight")
```

## Data Preparation

`massSight` works with the output of LC-MS experiments, which should
contain columns corresponding to:

1.  Compound ID
2.  Retention Time
3.  Mass to Charge Ratio
4.  (Optional) Average Intensity across all samples
5.  (Optional) Metabolite Name

``` r
data(hp1)
head(hp1) |>
  knitr::kable()
```

### The `massSight` Object

`massSight` creates and uses the `MSObject` class to store data and
results pertaining to individual LC-MS experiments. Prior to alignment,
LC-MS data frames or tibbles should be converted into an `MSObject`
using `create_ms_obj`:

``` r
data(hp1)
data(hp2)

ms1 <-
  create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity"
  )

ms2 <-
  create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity"
  )
```

An `MSObject` provides the following functions:

- `raw_df()` to access the experiment’s raw LC-MS data
- `isolated()` to access the experiment’s isolated metabolites, which is
  important for downstream alignment tasks
- `scaled_df()` to access the experiment’s scaled LC-MS data
- `consolidated()` to access the experiment’s consolidated data
- `metadata()` to access the experiment’s metadata

``` r
ms2 |>
  raw_df() |>
  head() |>
  knitr::kable(format = "simple")
```

| Compound_ID | Metabolite      |       RT |       MZ | Intensity |
|:------------|:----------------|---------:|---------:|----------:|
| cmp.3837    | C10 carnitine   | 7.261300 | 316.2479 | 638168.92 |
| cmp.3903    | C10:2 carnitine | 7.395033 | 312.2165 |  50418.96 |
| cmp.3749    | C12 carnitine   | 7.074067 | 344.2792 | 203210.69 |
| cmp.3756    | C12:1 carnitine | 7.105283 | 342.2635 | 363021.48 |
| cmp.3682    | C14 carnitine   | 6.926967 | 372.3107 |  93491.07 |
| cmp.3705    | C14:2 carnitine | 6.993833 | 368.2792 | 235545.00 |

## Alignment

### `auto_combine()`

Alignment is performed using `auto_combine()`

``` r
aligned <- auto_combine(
  ms1 = ms1,
  ms2 = ms2,
  iso_method = "dbscan"
)
```

More information on the `auto_combine()` function can be found in the
[package
documentation](https://omicseye.github.io/massSight/reference/auto_combine.html)

### `ml_match()`

``` r
ml_match_aligned <- ml_match(ms1, 
                             ms2, 
                             mz_thresh = 15, 
                             rt_thresh = 0.5, 
                             seed = 72)
```

### Plotting results from alignment

``` r
final_plots(aligned)
```

![](man/figures/final_plot_out.png)

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`
