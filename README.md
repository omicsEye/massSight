
<!-- README.md is generated from README.Rmd. Please edit that file -->

# massSight <img src="man/figures/massSight.png" align="right" height="300"/></a>

`massSight` is an R package for the alignment and scaling of LC-MS
metabolomics data.

- Citation: if you use `massSight`, please cite our manuscript: Chiraag
  Gohel and Ali Rahnavard. (2022). massSight: Metabolomics meta-analysis
  through multi-study data scaling, integration, and harmonization.
  <https://github.com/omicsEye/massSight>

## Contents

- [Dev Instructions](#dev-instructions)
- [Description](#description)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Input Files](#input-files-format)
- [Run](#run)
  - [Run a Demo](#run-a-demo)
- [Visualization](#visualization)
- [Example of real world
  applications](#example-of-real-world-applications)

## Dev Instructions

### Installation

1.  Clone/pull `massSight`
2.  Open the R project `massSight.Rproj`
3.  Build package using `devtools::build()`
4.  Install package using `devtools::install()`

### Usage/Testing

For now, manually add tsv/csv file and load data object into the global
environment using `read_csv` or `read_table`. Then run `auto_align` on
the created data objects.

For example:

``` r
library(massSight)

df1 <- readr::read_csv("HP-1.csv")
df2 <- readr::read_csv("HP-2.csv")
aligned_df <- auto_align(df1, df2)
```

## Description

## Requirements

## Installation

    devtools::install_github("omicsEye/massSight")

### Input files format

## Run

### Usage

Alignment is performed using `auto_align()`, and scaling is performed
using `auto_scale()`.

### Output Files

### Run a Demo

## Visualization

## Example of Real world applications
