
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
- [Requirements](#requirements)
- [Installation](#installation)
  - [Input Files](#input-files-format)
- [Run](#run)
  - [Run a Demo](#run-a-demo)
- [Visualization](#visualization)
- [Example of real world
  applications](#example-of-real-world-applications)
- [Dev Instructions](#dev-instructions)

## Description

## Requirements

## Installation

    devtools::install_github("omicsEye/massSight")

### Usage/Testing

Alignment is performed using `auto_align()`

``` r
library(massSight)

# example 1 inputs (small input for test)
View (df1_small)
View(df2_small)
aligned_df <- auto_align(df1_small, df2_small)


# example 2 inputs (small input for test)
aligned_df <- auto_align(df1, df2)
```

### Input files format

## Run

For now, manually add tsv/csv file and load data object into the global
environment using `read_csv` or `read_table`. Then run `auto_align` on
the created data objects.

``` r
# read inputs as data frame in R
user_df1 <- read.delim(
  "/path-to-file/profile1.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

user_df1 <- read.delim(
  "/path-to-file/profile2.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# run alignment function to combine two datasets
aligned_df <- auto_align(user_df1, user_df2)
```

### Output Files

## Visualization

``` r
final_plots(
  aligned_df$results_df_complete,
  aligned_df$smooth_for_plot,
  aligned_df$adjusted_df
)
```

## Example of Real world applications

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`
