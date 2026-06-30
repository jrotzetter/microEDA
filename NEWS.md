# microEDA 1.0.1

## Improvements

* `show_filter_history()` and `filter_features()` now use the term "features" instead of "taxa". This change more accurately reflects that rows may represent ASVs or OTUs that lack formal taxonomic classification.

# microEDA 1.0.0

* First official stable release on GitHub.
* Package is feature-complete for the initial release cycle.

## Features

### Core Functionality
- **`microEDA-class`**: Extends the `phyloseq-class` to store additional metadata and analysis history.
- **`microEDA-class` summarisation**: `show()`, `summary()`, and `show_filter_history()` print overviews of `microEDA` objects to the console.
- **MetaPhlAn support**: Load and merge MetaPhlAn profiles using `load_metaphlan()` and `join_mpa_profiles()`.
- **Feature filtering**: `filter_features()` filters taxa by abundance and prevalence, supporting both global and stratified (group-wise) operations.
- **Taxonomic agglomeration**: `agglomerate_taxa()` aggregates taxa at specified ranks; a vectorized, faster alternative to `phyloseq::tax_glom()`.
- **Taxonomy table utilities**: Functions to trim or add taxonomic prefixes (e.g., `k__`, `p__`) in a `taxonomyTable`.
- **Conversion to phyloseq**: `to_phyloseq()` converts `metaphlanProfile` and `microEDA` objects into standard `phyloseq` objects.
- **Taxonomic consistency check**: `check_taxonomic_consistency()` identifies inconsistencies within a `taxonomyTable`.
- **Presence lists**: `get_presence_list()` returns unique taxa present in each group defined by a sample metadata variable.
- **Taxa overlap analysis**: `get_taxa_overlaps()` computes overlaps and unique sets of taxa across sample groups, serving as a programmatic alternative to UpSet plots.

### Visualization Functions
`microEDA` provides a collection of `ggplot2`-based plotting functions for exploratory data analysis:

- **`plot_taxa_barchart()`**: Visualizes taxonomic abundance across samples, supporting aggregation at any rank and grouping by metadata.
- **`plot_taxa_heatmap()`**: Displays mean abundance and prevalence (detection frequency) of taxa within groups.
- **`plot_taxa_upset()`**: Visualizes shared and unique taxa across multiple sample groups (complemented by `get_taxa_overlaps()`).
- **`plot_taxa_sankey()`**: Illustrates hierarchical taxonomic relationships and abundance flow from higher to lower ranks.   
