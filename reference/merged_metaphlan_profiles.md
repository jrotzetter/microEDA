# Merged MetaPhlAn4 profiles data

Example data of multiple MetaPhlAn4 profiles that were merged into one
file in `metaphlanProfile` format.

## Usage

``` r
merged_metaphlan_profiles
```

## Format

### `merged_abundance_profiles`

An object of class
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
with slots `mpa_version` and `mpa_taxProfile`:

- mpa_version:

  A `Character` vector indicating the MetaPhlAn database version used.

- mpa_taxProfile:

  A data frame with 59 rows and 7 columns:

  `clade_name`: The taxonomic lineage of the taxon reported on this row.
  Taxon names are prefixed with one-letter codes to help indicate their
  rank.

  `SRS0144...`: The sample names for which sample-specific abundance
  profiles are found in a column.

## Source

<https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0>
