# microEDA

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
- **`microEDA-class` summarisation**: `show()`,
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`show_filter_history()`](https://jrotzetter.github.io/microEDA/reference/show_filter_history.md)
  print overviews of `microEDA` objects.
- **MetaPhlAn support**: Load and merge MetaPhlAn profiles using
  [`load_metaphlan()`](https://jrotzetter.github.io/microEDA/reference/load_metaphlan.md)
  and
  [`join_mpa_profiles()`](https://jrotzetter.github.io/microEDA/reference/join_mpa_profiles.md).
- **Feature filtering**:
  [`filter_features()`](https://jrotzetter.github.io/microEDA/reference/filter_features.md)
  allows filtering by abundance and prevalence, either globally or
  within groups defined by a metadata variable (stratified filtering).
- **Taxonomic agglomeration**:
  [`agglomerate_taxa()`](https://jrotzetter.github.io/microEDA/reference/agglomerate_taxa.md)
  aggregates taxa at specified taxonomic ranks. This implementation is
  significantly faster than
  [`phyloseq::tax_glom()`](https://rdrr.io/pkg/phyloseq/man/tax_glom.html).
- **Taxonomy table utilities**: Functions to trim or add taxonomic
  prefixes (e.g., `k__`, `p__`) in a `taxonomyTable`.
- **Conversion to phyloseq**:
  [`to_phyloseq()`](https://jrotzetter.github.io/microEDA/reference/to_phyloseq.md)
  converts `metaphlanProfile` and `microEDA` objects into standard
  `phyloseq` objects.
- **Taxonomic consistency check**:
  [`check_taxonomic_consistency()`](https://jrotzetter.github.io/microEDA/reference/check_taxonomic_consistency.md)
  identifies inconsistencies within a `taxonomyTable`.
- **Presence lists**:
  [`get_presence_list()`](https://jrotzetter.github.io/microEDA/reference/get_presence_list.md)
  returns unique taxa present in each group defined by a sample metadata
  variable.
- **Taxa overlap analysis**:
  [`get_taxa_overlaps()`](https://jrotzetter.github.io/microEDA/reference/get_taxa_overlaps.md)
  computes overlaps and unique sets of taxa across sample groups,
  serving as an alternative to UpSet plots.

### Visualization Functions

`microEDA` provides a collection of plotting functions, built on
ggplot2, and designed for exploratory data analysis of microbiome data.

#### 1. Taxonomic Composition Barplot

[`plot_taxa_barchart()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_barchart.md)
visualizes abundance of taxa across samples.

- **Aggregates** data at any taxonomic rank (e.g., Phylum, Genus).
- Optionally **Groups** samples by metadata variables (e.g., treatment,
  disease state).

#### 2. Mean Abundance & Prevalence Heatmap

[`plot_taxa_heatmap()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_heatmap.md)
displays a dual-metric view where each cell represents:

- **Mean abundance**: Average taxon abundance within a group.
- **Prevalence**: Proportion of samples in the group where the taxon is
  detected.

This dual metric helps distinguish consistently abundant taxa from those
that are sporadically present.

#### 3. Taxonomic Intersection UpSet Plot

[`plot_taxa_upset()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_upset.md)
visualizes shared and unique taxa across multiple sample groups.

- Complemented by
  [`get_taxa_overlaps()`](https://jrotzetter.github.io/microEDA/reference/get_taxa_overlaps.md)
  for programmatic access to intersection data.

#### 4. Taxonomic Flow Sankey Plot

[`plot_taxa_sankey()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_sankey.md)
illustrates hierarchical taxonomic relationships from higher (e.g.,
Phylum) to lower (e.g., Species) taxonomic ranks.

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

For more details please see
[`vignette("microEDA")`](https://jrotzetter.github.io/microEDA/articles/microEDA.md)
or the help pages in the documentation. Both are also available online
at <https://jrotzetter.github.io/microEDA/>.

## Future Directions

Potential areas for future development include:

- **Improved Color Scalability:** Exploring custom color palettes to
  better distinguish many groups.
- **Quality Control Extensions:** Investigating functions for data
  quality control, such as:
  - Contaminant removal and sequencing depth analysis.
  - Relative-abundance-aware metrics (e.g., zero-inflation rates).
  - Integration of outlier detection and rarefaction curves.

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
