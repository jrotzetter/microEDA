# Get or set the MetaPhlAn database version.

Extracts or updates the MetaPhlAn database version metadata stored in a
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
object. This version string reflects the reference database used during
profiling.

## Usage

``` r
mpa_version(object)

mpa_version(object) <- value
```

## Arguments

- object:

  An object of class `metaphlanProfile`.

- value:

  A `character` string with a new MetaPhlAn database version (used only
  in assignment).

## Value

The recorded MetaPhlAn database version (when getting), or updated
`metaphlanProfile` (when setting).

## See also

[`mpa_version,microEDA-method`](https://jrotzetter.github.io/microEDA/reference/info-accessors.md)
