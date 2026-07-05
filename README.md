
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
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

## Overview

`microEDA` is an R package for **exploratory data analysis of microbiome
data**, designed to answer the fundamental question: *“Who is in my
samples?”*

In addition to core visualization functions, the package provides some
utility functions to detect common data issues.

## Installation

The development version of microEDA can be installed from
[GitHub](https://github.com/jrotzetter/microEDA).

**Tip:** The [Introduction to
microEDA](https://jrotzetter.github.io/microEDA/articles/microEDA.html)
vignette is available online and can be viewed without installing the
package.

### Option 1: Using `pak` (Recommended)

This method is fastest and handles Bioconductor dependencies (like
`phyloseq`) automatically. Note that **vignettes are not included** with
this method.

``` r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("jrotzetter/microEDA")
```

### Option 2: Using `remotes` (Include Vignettes)

Select this method to include package vignettes for offline viewing.

***Important***: *Because `microEDA` depends on `phyloseq` (hosted on
Bioconductor), please install `phyloseq` first to avoid dependency
errors when using `remotes`.*

``` r
# 1. Install Bioconductor dependency
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

# 2. Install microEDA with vignettes
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
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

### Visualization Functions

`microEDA` provides a collection of plotting functions, built on
ggplot2, and designed for exploratory data analysis of microbiome data.

#### 1. Taxonomic Composition Barplot

`plot_taxa_barchart()` visualizes abundance of taxa across samples.

- **Aggregates** data at any taxonomic rank (e.g., Phylum, Genus).
- Optionally **Groups** samples by metadata variables (e.g., treatment,
  disease state).

#### 2. Mean Abundance & Prevalence Heatmap

`plot_taxa_heatmap()` displays a dual-metric view where each cell
represents:

- **Mean abundance**: Average taxon abundance within a group.
- **Prevalence**: Proportion of samples in the group where the taxon is
  detected.

This dual metric helps distinguish consistently abundant taxa from those
that are sporadically present.

#### 3. Taxonomic Intersection UpSet Plot

`plot_taxa_upset()` visualizes shared and unique taxa across multiple
sample groups.

- Complemented by `get_taxa_overlaps()` for programmatic access to
  intersection data.

#### 4. Taxonomic Flow Sankey Plot

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

For more details please see `vignette("microEDA")` or the help pages in
the documentation. Both are also available online at
<https://jrotzetter.github.io/microEDA/>.

## Planned features

- Custom color palette to increase available distinct colors is planned

## Citation

`microEDA` is a free, open-source project licensed under GPLv3. If you
find these functions useful for your publications, presentations, or
commercial workflows, please consider citing the package. This helps
track the project’s impact and ensures proper credit for the development
of these exploratory data analysis utilities.

**BibTeX:**

``` bibtex
@Manual{microEDA,
  title = {{microEDA}: Exploratory Microbiome Data Analysis and Visualization},
  author = {Jérémy Rotzetter},
  year = {2026},
  note = {R package version 1.0.1},
  url = {https://github.com/jrotzetter/microEDA},
}   
```

## Acknowledgments

`microEDA` is built on top of the
[phyloseq](https://github.com/joey711/phyloseq) infrastructure and
leverages the visualization capabilities of
[ggplot2](https://ggplot2.tidyverse.org/),
[ComplexUpset](https://github.com/krassowski/complex-upset), and
[ggsankeyfier](https://github.com/pepijn-devries/ggsankeyfier/).
Development of this package was made possible by the foundational work
of the authors and maintainers of these packages.

*Note: If you use `microEDA` in your research, please also consider
citing the underlying packages you utilize, particularly phyloseq.*
