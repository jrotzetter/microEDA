# Join multiple MetaPhlAn profiles

Combines multiple
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
objects into a single profile by joining their taxonomic abundance
tables on the `clade_name` column.

## Usage

``` r
join_mpa_profiles(...)
```

## Arguments

- ...:

  Multiple objects of class `metaphlanProfile` to merge.

## Value

An object of class `metaphlanProfile` containing the merged taxonomic
profiles. All input profiles must use the same MetaPhlAn database
version.

## Details

- Non-`metaphlanProfile` arguments are silently dropped with a warning.

- Profiles from different MetaPhlAn database versions will cause an
  error.

- Numeric suffixes will automatically be added to duplicate samples
  (e.g. Sample_1 -\> Sample_1.1, Sample_1.2).

- Missing values (`NA`) after joining are replaced with 0.0.

## Examples

``` r
if (FALSE) { # \dontrun{
mtp1 <- load_metaphlan("path/to/profile1.txt")
mtp2 <- load_metaphlan("path/to/profile2.txt")
merged_profile <- join_mpa_profiles(mtp1, mtp2)
} # }
```
