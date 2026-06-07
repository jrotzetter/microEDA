
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
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

## Overview

`microEDA` is an R package for **exploratory data analysis of microbiome
data**, designed to answer the fundamental question: *“Who is in my
samples?”*

In addition to core visualization tools, the package provides some
utility functions to detect common data issues.

## Installation

The development version of microEDA can be installed from
[GitHub](https://github.com/jrotzetter/microEDA) with:

``` r
# install.packages("pak")
pak::pak("jrotzetter/microEDA")
```

Or the `remotes` package:

``` r
# install.packages("remotes")
remotes::install_github("jrotzetter/microEDA", build_vignettes = TRUE)
```

## Currently Implemented Features

- **`microEDA-class`**: Extends the `phyloseq-class` to hold additional
  information.
- **`microEDA-class` summarisation**: `show()`, `summary()` and
  `show_filter_history()` print overviews of `microEDA` objects.
- **MetaPhlAn support**: Load and merge MetaPhlAn profiles using
  `load_metaphlan()` and `join_mpa_profiles()`.
- **Feature filtering**: `filter_features()` allows filtering by
  abundance and prevalence, either globally or within groups defined by
  a metadata variable (stratified filtering).
- **Taxonomic agglomeration**: `agglomerate_taxa()` aggregates taxa at
  specified taxonomic ranks. This implementation is significantly faster
  than `phyloseq::tax_glom()`.
- **Taxonomy table utilities**: Functions to trim or add taxonomic
  prefixes (e.g., `k__`, `p__`) in a `taxonomyTable`.
- **Conversion to phyloseq**: `to_phyloseq()` converts
  `metaphlanProfile` and `microEDA` objects into standard `phyloseq`
  objects.
- **Taxonomic consistency check**: `check_taxonomic_consistency()`
  identifies inconsistencies within a `taxonomyTable`.
- **Presence lists**: `get_presence_list()` returns unique taxa present
  in each group defined by a sample metadata variable.
- **Taxa overlap analysis**: `get_taxa_overlaps()` computes overlaps and
  unique sets of taxa across sample groups, serving as an alternative to
  UpSet plots.

## Visualization Functions

`microEDA` provides a collection of plotting functions, built on
ggplot2, and designed for exploratory data analysis of microbiome data.

### 1. Taxonomic Composition Barplot

`plot_taxa_barchart()` visualizes abundance of taxa across samples.

- **Aggregates** data at any taxonomic rank (e.g., Phylum, Genus).
- Optionally **Groups** samples by metadata variables (e.g., treatment,
  disease state).

### 2. Mean Abundance & Prevalence Heatmap

`plot_taxa_heatmap()` displays a dual-metric view where each cell
represents:

- **Mean abundance**: Average taxon abundance within a group.
- **Prevalence**: Proportion of samples in the group where the taxon is
  detected.

This dual metric helps distinguish consistently abundant taxa from those
that are sporadically present.

### 3. Taxonomic Intersection UpSet Plot

`plot_taxa_upset()` visualizes shared and unique taxa across multiple
sample groups.

- Complemented by `get_taxa_overlaps()` for programmatic access to
  intersection data.

### 4. Taxonomic Flow Sankey Plot

`plot_taxa_sankey()` illustrates hierarchical taxonomic relationships
from higher (e.g., Phylum) to lower (e.g., Species) taxonomic ranks.

- Supports display of abundance flow for a **single sample** or **mean
  abundance** across a group of samples.

## Example Usage

All functions are compatible with `phyloseq` objects, though the package
internal `microEDA-class` provides some additional information and
functionalities.

A `microEDA` object can be created by calling:

``` r
library(microEDA)
data(GlobalPatterns, package = "phyloseq")

me <- microEDA(GlobalPatterns)
```

In turn, a `microEDA` or `metaphlanProfile` object can be converted to
`phyloseq` with:

``` r
ps <- to_phyloseq(me)

identical(GlobalPatterns, ps)
#> [1] TRUE
```

## Planned features

- Custom color palette to increase available distinct colors is planned
