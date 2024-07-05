# massSight


- [Examples](#examples)
- [Description](#description)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
  - [The `massSight` Object](#the-masssight-object)
- [Alignment](#alignment)
  - [`auto_combine()`](#auto_combine)
  - [`ml_match()`](#ml_match)
- [Results](#results)
  - [Plotting results from alignment](#plotting-results-from-alignment)
  - [Using `massSight` to annotate unknown
    metabolites](#using-masssight-to-annotate-unknown-metabolites)
- [Dev Instructions](#dev-instructions)
  - [Installation](#installation-1)

<!-- README.md is generated from README.qmd. Please edit that file -->

<img src="man/figures/massSight.png" align="right" width="30%"/></a>

[![](https://zenodo.org/badge/608216683.svg)](https://zenodo.org/badge/latestdoi/608216683)

`massSight` is an R package for combining and scaling LC-MS metabolomics
data.

- Citation: if you use `massSight`, please cite our manuscript: Chiraag
  Gohel and Ali Rahnavard. (2023). massSight: Metabolomics meta-analysis
  through multi-study data scaling, integration, and harmonization.
  <https://github.com/omicsEye/massSight>

## Examples

Examples and extensive documentation can be found
[here](omicseye.github.io/massSight/)

## Description

## Installation

First, if you don’t have it installed, install `devtools` using:

``` r
install.packages("devtools")
```

Then, in an `R` console, run:

``` r
devtools::install_github("omicsEye/massSight")
```

You can then load the library using:

``` r
library(massSight)
```

## Data Preparation

`massSight` works with the output of LC-MS experiments, which should
contain columns corresponding to:

1.  Compound ID
2.  Retention Time
3.  Mass to Charge Ratio
4.  (Optional) Average Intensity across all samples
5.  (Optional) Metabolite Name

| Compound_ID      |       MZ |   RT | Intensity | Metabolite             |
|:-----------------|---------:|-----:|----------:|:-----------------------|
| 1.69_121.1014m/z | 121.1014 | 1.69 |  40329.32 | 1.2.4-trimethylbenzene |
| 3.57_197.0669m/z | 197.0669 | 3.57 | 117400.93 | 1,7-dimethyluric acid  |
| 7.74_282.1194m/z | 282.1194 | 7.74 |  16491.00 | 1-methyladenosine      |
| 5.27_166.0723m/z | 166.0723 | 5.27 |  22801.91 | 1-methylguanine        |
| 5.12_298.1143m/z | 298.1143 | 5.12 |  41602.96 | 1-methylguanosine      |
| 9.58_126.1028m/z | 126.1028 | 9.58 |   3004.32 | 1-methylhistamine      |

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
  ms1,
  ms2,
  rt_lower = -.5,
  rt_upper = .5,
  mz_lower = -15,
  mz_upper = 15,
  smooth_method = "gam",
  log = NULL
)
```

More information on the `auto_combine()` function can be found in the
[package
documentation](https://omicseye.github.io/massSight/reference/auto_combine.html)

### `ml_match()`

The `ml_match()` function is an alternative method for merging LC-MS
experiments with semi-annotated data sets.

``` r
ml_match_aligned <- ml_match(
  ms1,
  ms2,
  mz_thresh = 15,
  rt_thresh = 0.5,
  prob_thresh = .5,
  seed = 72
)
```

## Results

Results from an alignment function are stored as a `MergedMSObject`.
This object contains the following slots:

- `all_matched()`: All of the final matched metabolites between the two
  datasets. This is the main result of the various matching functions.

``` r
all_matched(aligned) |>
  head() |>
  knitr::kable()
```

| rep_Compound_ID | rep_RT | rep_MZ | rep_Intensity | rep_Metabolite | Compound_ID_hp1 | Compound_ID_hp2 | Metabolite_hp1 | Metabolite_hp2 | RT_hp1 | RT_hp2 | MZ_hp1 | MZ_hp2 | Intensity_hp1 | Intensity_hp2 |
|:---|---:|---:|---:|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|
| 0.63_74.0995m/z | 0.630000 | 74.0995 | 0.40 |  | 0.63_74.0995m/z | cmp.1 |  |  | 0.63 | 1.003550 | 74.0995 | 74.09722 | 0.40 | 5683858.19 |
| cmp.10 | 1.003550 | 429.4048 | 14695.70 |  | NA | cmp.10 | NA |  | NA | 1.003550 | NA | 429.40477 | NA | 14695.70 |
| 1.14_262.0378m/z | 1.140000 | 262.0378 | 17383.20 |  | 1.14_262.0378m/z | cmp.100 |  |  | 1.14 | 1.199700 | 262.0378 | 262.03783 | 17383.20 | 20343.60 |
| cmp.1000 | 1.918017 | 540.5345 | 38459.42 |  | NA | cmp.1000 | NA |  | NA | 1.918017 | NA | 540.53450 | NA | 38459.42 |
| 1.74_348.0645m/z | 1.740000 | 348.0645 | 59393.19 |  | 1.74_348.0645m/z | cmp.1001 |  |  | 1.74 | 1.918017 | 348.0645 | 348.06456 | 59393.19 | 41235.47 |
| 1.86_478.1970m/z | 1.860000 | 478.1970 | 132193.24 |  | 1.86_478.1970m/z | cmp.1002 |  |  | 1.86 | 1.918017 | 478.1970 | 478.19718 | 132193.24 | 35052.10 |

- `iso_matched()`: The matched isolated metabolites between the two
  datasets.

``` r
iso_matched(aligned) |>
  head() |>
  knitr::kable()
```

| df1 | RT | MZ | Intensity | df2 | RT_2 | MZ_2 | Intensity_2 | delta_RT | smooth_rt | srt | delta_MZ | smooth_mz | smz | sintensity |
|:---|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7.98_76.0402m/z | 7.98 | 76.0402 | 13067.16 | cmp.4356 | 8.193933 | 76.04010 | 27158.674 | 0.2139333 | 0.1462600 | 7.848138 | -0.0001027 | -9.2e-06 | 76.04021 | 7832.970 |
| 8.14_77.0799m/z | 8.14 | 77.0799 | 9291.18 | cmp.4388 | 8.283100 | 77.07982 | 11295.184 | 0.1431000 | 0.1617242 | 7.992453 | -0.0000786 | -8.9e-06 | 77.07991 | 3625.042 |
| 4.79_79.0220m/z | 4.79 | 79.0220 | 2356.37 | cmp.2649 | 4.910117 | 79.02198 | 6818.932 | 0.1201167 | 0.1440289 | 4.627987 | -0.0000210 | -8.1e-06 | 79.02201 | 2327.172 |
| 7.50_80.1318m/z | 7.50 | 80.1318 | 84942.05 | cmp.4025 | 7.635767 | 80.13173 | 11765.610 | 0.1357667 | 0.0897390 | 7.419991 | -0.0000662 | -7.7e-06 | 80.13181 | 3757.301 |
| 7.86_84.9120m/z | 7.86 | 84.9120 | 10189.64 | cmp.4247 | 8.015617 | 84.91193 | 9442.532 | 0.1556167 | 0.1332109 | 7.740726 | -0.0000657 | -6.0e-06 | 84.91201 | 3097.303 |
| 8.75_84.9605m/z | 8.75 | 84.9605 | 240071.83 | cmp.4584 | 8.978550 | 84.96037 | 287055.188 | 0.2285500 | 0.1999628 | 8.559125 | -0.0001285 | -5.9e-06 | 84.96051 | 62125.356 |

### Plotting results from alignment

The `final_plots()` function returns plots containing information on RT
and MZ drift for pre isolation, isolation, and final matching results.
These plots can be used for diagnostic purposes.

``` r
plots <- final_plots(aligned,
  rt_lim = c(-.5, .5),
  mz_lim = c(-15, 15)
)
plots
```

![](man/figures/final_plot_out.png)

This plot can be saved locally using `ggsave()` from the `ggplot2`
package:

``` r
ggplot2::ggsave(
  filename = "plot.png",
  plot = plots
)
```

### Using `massSight` to annotate unknown metabolites

``` r
merged_df <- all_matched(aligned)
hp2_annotated <- merged_df |>
  dplyr::select(Compound_ID_hp2, rep_Metabolite) |>
  dplyr::inner_join(hp2, by = c("Compound_ID_hp2" = "Compound_ID"))

hp2_annotated |> 
  dplyr::filter(rep_Metabolite != "") |>
  dplyr::arrange(rep_Metabolite) |>
  head(10) |>
  knitr::kable()
```

| Compound_ID_hp2 | rep_Metabolite | Metabolite | RT | MZ | Intensity |
|:---|:---|:---|---:|---:|---:|
| cmp.2157 | 1,7-dimethyluric acid |  | 3.817950 | 197.0669 | 140175.90 |
| cmp.4168 | 1-methyladenosine |  | 7.864050 | 282.1194 | 43167.05 |
| cmp.2810 | 1-methylguanine |  | 5.361300 | 166.0723 | 25898.35 |
| cmp.2782 | 1-methylguanosine |  | 5.267700 | 298.1145 | 70138.72 |
| cmp.785 | 1.2.4-trimethylbenzene |  | 1.805967 | 121.1014 | 13453.95 |
| cmp.2750 | 3-(N-acetyl-L-cystein-S-yl) acetaminophen |  | 5.124100 | 313.0850 | 38069.50 |
| cmp.4740 | 3-methylhistidine |  | 10.289133 | 170.0923 | 490315.34 |
| cmp.2091 | 4-acetamidobutanoate |  | 3.595050 | 146.0811 | 78624.40 |
| cmp.2828 | 5-acetylamino-6-amino-3-methyluracil |  | 5.424667 | 199.0824 | 29391.98 |
| cmp.2446 | 6.8-dihydroxypurine |  | 4.544583 | 153.0407 | 10692.24 |

Here, `rep_Metabolite` is the metabolite name from the reference
dataset.

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`
