# A S4 class extending phyloseq for microbiome data analysis

This class inherits from
[`phyloseq`](https://rdrr.io/pkg/phyloseq/man/phyloseq.html) and adds a
slot for additional information, such as lowest taxonomic rank,
MetaPhlAn database version, applied transformations and/or filters.

`show()` displays a summary of a `microEDA` object, including the
standard `phyloseq` output and information from the `info` slot.

## Usage

``` r
microEDA(tax_profile, metadata = NULL, sample_column = NULL)

# S4 method for class 'metaphlanProfile'
microEDA(tax_profile, metadata = NULL, sample_column = NULL)

# S4 method for class 'phyloseq'
microEDA(tax_profile, metadata = NULL, sample_column = NULL)

# S4 method for class 'microEDA'
show(object)
```

## Arguments

- tax_profile:

  An object of class `metaphlanProfile` or `phyloseq`.

- metadata:

  Optional `data.frame` with sample metadata.

- sample_column:

  `Character` string specifying the column name in `metadata` that
  contains sample identifiers (e.g., "Samples").

- object:

  A `microEDA` object.

## Value

An object of class `microEDA`.

## Details

When passed a
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md),
the profile will be returned at the lowest available taxonomic rank,
which usually is at strain rank. To agglomerate the profile at a higher
rank, please see
[`agglomerate_taxa()`](https://jrotzetter.github.io/microEDA/reference/agglomerate_taxa.md).

## Slots

- `info`:

  A `list` containing metadata.

- `otu_table`:

  A `matrix` of class `otu_table` with taxa as rows and samples as
  columns.

- `tax_table`:

  A named `character matrix` of class `taxonomyTable`, with taxon IDs as
  row names and taxonomic ranks as columns.

- `sam_data`:

  An object of class `sample_data` containing sample metadata.

- `phy_tree`:

  An object of class `phylo` representing a phylogenetic tree.

- `refseq`:

  An object of a class inheriting from `XStringSet-class` (e.g.,
  `DNAStringSet`) containing biological sequences. See
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html).

## Examples

``` r
data("GlobalPatterns", package = "phyloseq")
me <- microEDA(GlobalPatterns)
```
