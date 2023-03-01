
<!-- README.md is generated from README.Rmd. Please edit that file -->

# massSight <img src="man/figures/massSight.png" align="right" height="300"/></a>

`massSight` is an R package for the alignment and scaling of LC-MS
metabolomics data.

- Citation: if you use `massSight`, please cite our manuscript: Chiraag
  Gohel and Ali Rahnavard. (2023). massSight: Metabolomics meta-analysis
  through multi-study data scaling, integration, and harmonization.
  <https://github.com/omicsEye/massSight>

## Contents

- [Description](#description)
- [Dev Instructions](#dev-instructions)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Input Files](#input-files-format)
- [Run](#run)
  - [Run a Demo](#run-a-demo)
- [Visualization](#visualization)
- [Example of real world
  applications](#example-of-real-world-applications)

## Description

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`

### Usage/Testing

Alignment is performed using `auto_align()`

``` r
library(massSight)

# example 1 inputs (small input for test)
View (df1_small)
View(df2_samll)
aligned_df <- auto_align(df1_small, df2_small)


# example 2 inputs (small input for test)
aligned_df <- auto_align(df1, df2)
```

## Requirements

## Installation

    devtools::install_github("omicsEye/massSight")

### Input files format

## Run

For now, manually add tsv/csv file and load data object into the global
environment using `read_csv` or `read_table`. Then run `auto_align` on
the created data objects.

``` r
# read inputs as dataframe in R
user_df1 <- readr::read_csv("/path-to-file/profile1.tsv")
user_df2 <- readr::read_csv("/path-to-file/profile2.tsv")


# run alignment function to combine two datasets
aligned_df <- auto_align(user_df1, user_df2)
```

### Output Files

## Visualization

``` r
final_plots(
  res$results_df_complete,
  res$smooth_for_plot,
  res$adjusted_df
)
```

## Example of Real world applications
