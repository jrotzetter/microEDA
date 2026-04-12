
<!-- README.md is generated from README.Rmd. Please edit that file -->

# microEDA

<!-- badges: start -->

[![GitHub
Release](https://img.shields.io/github/release/jrotzetter/microEDA?include_prereleases=&sort=semver&color=blue)](https://github.com/jrotzetter/microEDA/releases/ "View releases")
[![License](https://img.shields.io/badge/License-GPLv3-blue)](#license "View license summary")
[![Issues -
microEDA](https://img.shields.io/github/issues/jrotzetter/microEDA)](https://github.com/jrotzetter/microEDA/issues "View open issues")
[![Made with
R](https://img.shields.io/badge/R-4.5.3-blue?logo=r&logoColor=white)](https://cran.r-project.org/ "Go to CRAN homepage")
[![Made with
R](https://img.shields.io/badge/RStudio-2026.01.1_Build_403-blue?logo=rstudio&logoColor=white)](https://posit.co/products/open-source/rstudio/ "Go to RSTUDIO IDE homepage")

<!-- badges: end -->

## Overview

R package for microbiome exploratory data analysis. Currently under
active development.

## Installation

You can install the development version of microEDA from
[GitHub](https://github.com/jrotzetter/microEDA) with:

``` r
# install.packages("devtools")
devtools::install_github("jrotzetter/microEDA", build_vignettes = TRUE)
```

Alternatively you can also use the `pak` package:

``` r
# install.packages("pak")
pak::pak("jrotzetter/microEDA")
```

Or the `remotes` package:

``` r
# install.packages("remotes")
remotes::install_github("jrotzetter/microEDA", build_vignettes = TRUE)
```

## Planned features

`microEDA` is planned to provide a collection of plotting functions,
built on ggplot2, and designed for exploratory analysis of microbiome
data.

1.  Taxonomic Composition Barplot
    - Visualize the relative abundance of taxa within samples
    - Aggregate data at any taxonomic rank (e.g., Phylum, Genus)
    - Group samples by metadata variables (e.g., treatment, disease
      state)
    - Custom color palette to increase available distinct colors is
      planned
2.  Mean Abundance & Prevalence Heatmap Table
    - Each cell will display:
      - Mean Relative Abundance of a taxon within a group
      - Prevalence indicating the proportion of samples within a group
        where a taxon is detected
    - This combined view helps distinguish between consistently abundant
      taxa and those that are highly variable
3.  Taxonomic Intersection UpSet Plot
    - Show the number of shared and unique taxa between any number of
      sample groups
    - Separate function to obtain list of intersections will be provided
4.  Taxonomic Flow Sankey Plot
    - Visualize the hierarchical relationship from higher (e.g., Phylum)
      to lower (e.g., Species) taxonomic ranks
    - Either plot the relative abundance “flow” for a single sample or
      the mean abundance across a group of samples
