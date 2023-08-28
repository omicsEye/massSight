# massSight

- [Examples](#examples)
- [Description](#description)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
  - [The `massSight` Object](#ms-obj)
- [Alignment](#align)
  - [`auto_combine()`](#auto_combine)
  - [`ml_match()`](#ml_match)
  - [Results](#results)
  - [Plotting results from alignment](#plotting-results-from-alignment)
- [Dev Instructions](#dev-instrutions)
  - [Installation](#installation-1)

<!-- README.md is generated from README.qmd. Please edit that file -->

<img src="man/figures/massSight.png" align="right" width="30%"/></a>

<div>

[![](https://zenodo.org/badge/608216683.svg)](https://zenodo.org/badge/latestdoi/608216683)

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

First, if you don’t have it installed, install `devtools` using:

``` r
install.packages("devtools")
```

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
  ms1 = ms1,
  ms2 = ms2
)
#> Numbers of matched/kept features: 370
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

### Results

Results from an alignment function are stored as a `MergedMSObject`.
This object contains the following slots:

- `all_matched()`: All of the final matched metabolites between the two
  datasets. This is the main result of the various matching functions.

``` r
all_matched(aligned) |>
  head() |>
  knitr::kable()
```

| df1              |   RT |       MZ |   Intensity | Metabolite              | df2      |     RT_2 |     MZ_2 |  Intensity_2 | Metabolite_2 |
|:-----------------|-----:|---------:|------------:|:------------------------|:---------|---------:|---------:|-------------:|:-------------|
| 6.91_267.1700m/z | 6.91 | 267.1700 |    39354.33 | atenolol                | cmp.3676 | 6.918050 | 267.1700 |     140703.1 |              |
| 7.29_703.5736m/z | 7.29 | 703.5736 | 34050464.62 | C16:0 SM                | cmp.3861 | 7.310333 | 703.5739 | 102344039\.7 |              |
| 7.68_480.3445m/z | 7.68 | 480.3445 |   130447.16 | C16:1 LPC plasmalogen   | cmp.4075 | 7.734783 | 480.3446 |     558532.5 |              |
| 7.71_524.3709m/z | 7.71 | 524.3709 | 10367615.65 | C18:0 LPC               | cmp.4096 | 7.770433 | 524.3710 |   24208681.0 |              |
| 7.22_731.6045m/z | 7.22 | 731.6045 |  8193010.61 | C18:0 SM                | cmp.3813 | 7.234550 | 731.6050 |   19793962.6 |              |
| 7.86_508.3760m/z | 7.86 | 508.3760 |   281006.87 | C18:1 LPC plasmalogen_A | cmp.4217 | 7.944300 | 508.3760 |     577530.9 |              |

- `iso_matched()`: The matched isolated metabolites between the two
  datasets.

``` r
iso_matched(aligned) |> 
  head() |>
  knitr::kable()
```

| df1              |   RT |       MZ | Intensity | df2      |     RT_2 |     MZ_2 | Intensity_2 |   delta_RT | smooth_rt |      srt |   delta_MZ | smooth_mz |      smz | sintensity |
|:-----------------|-----:|---------:|----------:|:---------|---------:|---------:|------------:|-----------:|----------:|---------:|-----------:|----------:|---------:|-----------:|
| 8.19_138.0525m/z | 8.19 | 138.0525 |  24972.37 | cmp.4390 | 8.287550 | 138.0525 |    56384.44 |  0.0975500 | 0.0839866 | 8.109608 |  0.0000451 | -4.95e-05 | 138.0525 |  17081.178 |
| 6.60_159.0916m/z | 6.60 | 159.0916 |  17528.13 | cmp.3556 | 6.699617 | 159.0916 |    11669.18 |  0.0996167 | 0.0219970 | 6.577464 |  0.0000081 | -1.05e-05 | 159.0916 |   4198.344 |
| 6.42_163.0991m/z | 6.42 | 163.0991 |   1382.40 | cmp.3431 | 6.445517 | 163.0991 |     9093.37 |  0.0255167 | 0.0274847 | 6.391444 | -0.0000376 | -3.40e-06 | 163.0991 |   3361.912 |
| 6.03_182.1903m/z | 6.03 | 182.1903 |  14936.08 | cmp.3244 | 6.057683 | 182.1903 |    19013.68 |  0.0276833 | 0.0419943 | 5.986692 | -0.0000203 |  2.86e-05 | 182.1903 |   6485.731 |
| 7.11_185.0766m/z | 7.11 | 185.0766 |   9829.32 | cmp.3753 | 7.096367 | 185.0766 |    42251.75 | -0.0136333 | 0.0213218 | 7.089485 | -0.0000241 |  3.32e-05 | 185.0766 |  13209.387 |
| 5.68_186.1318m/z | 5.68 | 186.1318 |  46177.82 | cmp.3016 | 5.660933 | 186.1317 |    61364.89 | -0.0190667 | 0.0528671 | 5.625650 | -0.0001496 |  3.48e-05 | 186.1318 |  18418.984 |

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
