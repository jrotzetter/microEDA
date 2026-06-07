# Load a MetaPhlAn taxonomic profile from a file

Reads a tab-separated MetaPhlAn output file and constructs a
[metaphlanProfile](https://jrotzetter.github.io/microEDA/reference/metaphlanProfile-class.md)
object. The input file should contain taxonomic abundance data with
columns for clade name and sample abundances.

## Usage

``` r
load_metaphlan(file_path)
```

## Arguments

- file_path:

  A `character` string containing the path to the file to read.

## Value

An object of class `metaphlanProfile`.

## Examples

``` r
file_path <- system.file("extdata", "merged_abundance_table.txt", package = "microEDA")

mtp <- load_metaphlan(file_path)
```
