# SRS014476-Supragingival_plaque data

Example data of a single MetaPhlAn4 output file in `metaphlanProfile`
format, containing the computed taxon abundances.

## Usage

``` r
single_metaphlan_profile
```

## Format

### `single_abundance_profile`

An object of class
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
with slots `mpa_version` and `mpa_taxProfile`:

- mpa_version:

  A `Character` vector indicating the MetaPhlAn database version used.

- mpa_taxProfile:

  A data frame with 20 rows and 2 columns:

  `clade_name`: The taxonomic lineage of the taxon reported on this row.
  Taxon names are prefixed with one-letter codes to help indicate their
  rank.

  `single_sample_profile`: Each taxon's relative abundance in %.

## Source

<https://github.com/biobakery/MetaPhlAn/releases/tag/4.1.0>
