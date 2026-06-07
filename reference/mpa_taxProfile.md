# Get or set the MetaPhlAn taxonomic profile

Extracts or updates the core taxonomic abundance data from a
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
object. The taxonomic profile is stored as a `data.frame` containing
taxonomic identifiers (clade_name) and relative abundances per sample.

## Usage

``` r
mpa_taxProfile(object)

mpa_taxProfile(object) <- value
```

## Arguments

- object:

  An object of class `metaphlanProfile`.

- value:

  A `data.frame` with a MetaPhlAn taxonomic profile. Taxonomic lineages
  should be stored in a 'clade_name' column.

## Value

The taxonomic profile (when getting), or updated `metaphlanProfile`
(when setting).
