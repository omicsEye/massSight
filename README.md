
<!-- README.md is generated from README.qmd. Please edit that file -->

# massSight

<img src="man/figures/massSight.png" align="right" width="30%"/></a>

`massSight` is an R package for the alignment and scaling of LC-MS
metabolomics data.

- Citation: if you use `massSight`, please cite our manuscript: Chiraag
  Gohel and Ali Rahnavard. (2023). massSight: Metabolomics meta-analysis
  through multi-study data scaling, integration, and harmonization.
  <https://github.com/omicsEye/massSight>

## Contents

- [Description](#description)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
  - [The `massSight` Object](#ms-obj)
- [Alignment](#align)
- [Visualization](#visualization)
- [Example of real world
  applications](#example-of-real-world-applications)
- [Dev Instructions](#dev-instructions)

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

| Compound_ID |       MZ |       RT | Intensity | Metabolite            |
|:------------|---------:|---------:|----------:|:----------------------|
| CMP-2758    | 197.0665 | 3.401200 |  74498.87 | 1,7-dimethyluric acid |
| CMP-4802    | 282.1189 | 8.627617 |  25684.75 | 1-methyladenosine     |
| CMP-3329    | 166.0720 | 6.069083 |  28585.37 | 1-methylguanine       |
| CMP-3294    | 298.1139 | 5.913067 |  61491.94 | 1-methylguanosine     |
| CMP-5077    | 137.0707 | 9.114433 | 112107.07 | 1-methylnicotinamide  |
| CMP-1830    | 347.2210 | 2.017650 |   1291.29 | 21-deoxycortisol      |

### The `massSight` Object

`massSight` creates and uses the `MSObject` class to store data and
results pertaining to individual LC-MS experiments. Prior to alignment,
LC-MS data frames or tibbles should be converted into an `MSObject`
using `create_ms_obj`:

``` r
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
  kableExtra::kbl(format = "simple")
```

| Compound_ID      |   RT |       MZ | Intensity |
|:-----------------|-----:|---------:|----------:|
| 1.68_121.1013m/z | 1.68 | 121.1013 | 50543.745 |
| 3.53_197.0667m/z | 3.53 | 197.0667 | 21948.556 |
| 7.81_282.1190m/z | 7.81 | 282.1190 | 14220.869 |
| 5.29_166.0721m/z | 5.29 | 166.0721 | 46434.807 |
| 5.16_298.1139m/z | 5.16 | 298.1139 | 60616.582 |
| 9.77_126.1026m/z | 9.77 | 126.1026 |  4435.973 |

## Alignment

Alignment is performed using `auto_align()`

``` r
aligned <- auto_align(ms1 = ms1, 
                      ms2 = ms2, 
                      iso_method = "dbscan")
```

More information on the `auto_align()` function can be found in the
[package
documentation](https://omicseye.github.io/massSight/reference/auto_align.html)

### Plotting results from alignment

``` r
final_plots(aligned)
```

![](man/figures/final_plot_out.png)

### Input files format

## Run

For now, manually add tsv/csv file and load data object into the global
environment using `read_csv` or `read_table`. Then run `auto_align` on
the created data objects.

``` r
# read inputs as data frame in R
user_df1 <- read.delim(
  "/path-to-file/profile1.tsv",
  sep = "\t",
  header = T,
  fill = F,
  comment.char = "",
  check.names = F,
  row.names = 1
)

user_df1 <- read.delim(
  "/path-to-file/profile2.tsv",
  sep = "\t",
  header = T,
  fill = F,
  comment.char = "",
  check.names = F,
  row.names = 1
)

# run alignment function to combine two datasets
aligned <- auto_align(user_df1, user_df2)
```

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`
