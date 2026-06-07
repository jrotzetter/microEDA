# Create an UpSet Plot of Shared and Unique Taxa

Generates a ggplot2-based UpSet plot to visualize the intersection
patterns of taxa across different sample groups. This provides a
scalable alternative to Venn diagrams for microbiome data.

## Usage

``` r
plot_taxa_upset(
  me,
  tax_rank,
  group_var,
  plot_title = NULL,
  show_names = FALSE,
  group_labels = NULL,
  min_abundance = 0,
  min_prevalence = 0,
  filter_by_group = FALSE,
  text_size = 2,
  ...
)
```

## Arguments

- me:

  A `microEDA` or `phyloseq` object containing OTU table, taxonomy, and
  sample data.

- tax_rank:

  `Character` string specifying the taxonomic rank to plot (e.g.,
  "Phylum", "Family"). Must be a valid rank in the taxonomy table.

- group_var:

  `Character` string indicating a sample metadata variable to group
  samples.

- plot_title:

  Optional `character` string for the plot title. If `NULL`, a default
  title is generated based on `tax_rank`.

- show_names:

  `Logical`. If `TRUE`, displays taxon names on intersection bars with
  italicized formatting. Default is `FALSE`.

- group_labels:

  Named character `vector` mapping old group names to new labels (e.g.,
  `c("Old" = "New")`).

- min_abundance:

  `Numeric` value. Minimum abundance threshold for a feature to be
  retained. Must be non-negative. Features with abundance below this are
  considered absent.

- min_prevalence:

  `Numeric` value. Minimum prevalence required for retention. If value
  is \< 1, interpreted as proportion of samples; otherwise, as absolute
  number of samples.

- filter_by_group:

  `Logical`. If `FALSE` (default), filtering is applied globally across
  all samples even if `group_var` is specified. This allows using
  `group_var` for stratification in plotting without affecting the
  filtering scope. If `TRUE`, filtering is applied within each group
  defined by `group_var` and `group_requirement`. See
  [filter_features](https://jrotzetter.github.io/microEDA/reference/filter_features.md)
  for more details on filtering arguments.

- text_size:

  `Numeric`. Font size for taxon labels when `show_names = TRUE`.
  Default is 2.

- ...:

  Additional arguments for fine-tuning. Can include:

  - `abundance_criterion`: `Character` string. Criterion to use for
    filtering:

    `prevalence`:

    :   Retain features present in at least `min_prevalence` samples
        (within group if `group_var` is used) and with abundance \>=
        `min_abundance` in those samples.

    `mean`:

    :   Also requires that the mean abundance across samples (or group)
        is \>= `min_abundance`.

    Default: `"prevalence"`.

  - `group_requirement`: `Character` string. When `group_var` is
    specified, determines whether the filter criterion must be met in
    `"any"` group or `"all"` groups. Default: `"any"`.

  - `keep_filtered`: Whether to keep filtered out taxa as "Other"
    (Default: `TRUE`). Takes effect only for `microEDA` objects.

  - `rm_missing`: `Logical.` If `TRUE`, removes taxa with
    missing/unclassified entries at the specified rank. If `FALSE`,
    fills missing values by propagating the last known ancestor,
    labeling them as "Unclassified *Last_Known_Parent_Clade*" (e.g.,
    "Unclassified Enterobacteriaceae"). (Default: `FALSE`).

## Value

A `ggplot` object.

## Details

The function performs the following steps:

1.  Validates the input `me` object.

2.  Filters taxa by `min_abundance`, `min_prevalence`.

3.  Aggregates counts by `tax_rank`.

4.  Determines taxon presence/absence per sample group.

5.  Constructs an UpSet plot using the `ComplexUpset` package.

Taxon names are formatted in italics when `show_names = TRUE`.

## Examples

``` r
data(GlobalPatterns, package = "phyloseq")
plot_taxa_upset(GlobalPatterns, "Phylum", "SampleType")
#> Warning: Both 'min_abundance' and 'min_prevalence' are 0. No filtering will be applied.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ComplexUpset package.
#>   Please report the issue at
#>   <https://github.com/krassowski/complex-upset/issues>.

```
