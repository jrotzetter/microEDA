# Changelog

## microEDA 1.0.1

### Improvements

- [`show_filter_history()`](https://jrotzetter.github.io/microEDA/reference/show_filter_history.md)
  and
  [`filter_features()`](https://jrotzetter.github.io/microEDA/reference/filter_features.md)
  now use the term “features” instead of “taxa”. This change more
  accurately reflects that rows may represent ASVs or OTUs that lack
  formal taxonomic classification.

## microEDA 1.0.0

- First official stable release on GitHub.
- Package is feature-complete for the initial release cycle.

### Features

#### Core Functionality

- **`microEDA-class`**: Extends the `phyloseq-class` to store additional
  metadata and analysis history.
- **`microEDA-class` summarisation**: `show()`,
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`show_filter_history()`](https://jrotzetter.github.io/microEDA/reference/show_filter_history.md)
  print overviews of `microEDA` objects to the console.
- **MetaPhlAn support**: Load and merge MetaPhlAn profiles using
  [`load_metaphlan()`](https://jrotzetter.github.io/microEDA/reference/load_metaphlan.md)
  and
  [`join_mpa_profiles()`](https://jrotzetter.github.io/microEDA/reference/join_mpa_profiles.md).
- **Feature filtering**:
  [`filter_features()`](https://jrotzetter.github.io/microEDA/reference/filter_features.md)
  filters taxa by abundance and prevalence, supporting both global and
  stratified (group-wise) operations.
- **Taxonomic agglomeration**:
  [`agglomerate_taxa()`](https://jrotzetter.github.io/microEDA/reference/agglomerate_taxa.md)
  aggregates taxa at specified ranks; a vectorized, faster alternative
  to
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
  serving as a programmatic alternative to UpSet plots.

#### Visualization Functions

`microEDA` provides a collection of `ggplot2`-based plotting functions
for exploratory data analysis:

- **[`plot_taxa_barchart()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_barchart.md)**:
  Visualizes taxonomic abundance across samples, supporting aggregation
  at any rank and grouping by metadata.
- **[`plot_taxa_heatmap()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_heatmap.md)**:
  Displays mean abundance and prevalence (detection frequency) of taxa
  within groups.
- **[`plot_taxa_upset()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_upset.md)**:
  Visualizes shared and unique taxa across multiple sample groups
  (complemented by
  [`get_taxa_overlaps()`](https://jrotzetter.github.io/microEDA/reference/get_taxa_overlaps.md)).
- **[`plot_taxa_sankey()`](https://jrotzetter.github.io/microEDA/reference/plot_taxa_sankey.md)**:
  Illustrates hierarchical taxonomic relationships and abundance flow
  from higher to lower ranks.
