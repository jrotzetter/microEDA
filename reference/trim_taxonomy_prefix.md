# Remove taxonomic prefixes from tax_table columns

Removes QIIME-style two-underscore prefixes (e.g., "k\_\_", "p\_\_")
from values in a taxonomy matrix or data.frame if present. The function
is idempotent: repeated application has no additional effect.

## Usage

``` r
trim_taxonomy_prefix(
  me,
  prefix_map = c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__", Family =
    "f__", Genus = "g__", Species = "s__", Strain = "t__")
)
```

## Arguments

- me:

  A `character` matrix, `data.frame`, `microEDA` or a `phyloseq`-class
  object containing a `tax_table`. If a `microEDA` or `phyloseq` object
  is passed, its taxonomyTable is modified.

- prefix_map:

  A named `character` vector mapping taxonomic rank names (e.g.,
  "Kingdom", "Phylum") to their corresponding prefixes (e.g., "k\_\_",
  "p\_\_"). Defaults to standard QIIME-style prefixes including "Strain"
  ("t\_\_").

## Value

Returns the input object with prefixes removed:

- For matrices/data.frames: returns a modified `character` matrix or
  data.frame.

- For `microEDA`/`phyloseq` objects: returns the updated
  `microEDA`/`phyloseq` object with cleaned `tax_table`.

## Details

Supports character matrices, data.frames, `microEDA` and `phyloseq`
objects. Only modifies cells in columns matching the keys of
`prefix_map` and only when the value starts with the corresponding
prefix. Missing (`NA`), empty, or unprefixed values are left unchanged.

- Only columns in `me` that match names in `prefix_map` are processed.

- Prefix removal is exact: e.g., "k\_\_" is only removed from "Kingdom"
  column if it starts the value.

- Uses efficient string processing (e.g., via `stringi`) to detect and
  remove prefixes.

- Preserves `NA`, `""`, and whitespace-only entries.

## See also

[`add_taxonomy_prefix()`](https://jrotzetter.github.io/microEDA/reference/add_taxonomy_prefix.md)
to add prefixes.

[`tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)
to access taxonomy data.

## Examples

``` r
# Example 1: Basic usage with data.frame
prefixed_tax_data <- data.frame(
  Kingdom = c("k__Bacteria", "k__Archaea"),
  Genus = c("g__Escherichia", "g__Methanobrevibacter"),
  Species = c("s__Escherichia coli", "s__Methanobrevibacter smithii"),
  stringsAsFactors = FALSE
)
prefixed_tax_data
#>       Kingdom                 Genus                       Species
#> 1 k__Bacteria        g__Escherichia           s__Escherichia coli
#> 2  k__Archaea g__Methanobrevibacter s__Methanobrevibacter smithii
trimmed <- trim_taxonomy_prefix(prefixed_tax_data)
trimmed
#>    Kingdom              Genus                    Species
#> 1 Bacteria        Escherichia           Escherichia coli
#> 2  Archaea Methanobrevibacter Methanobrevibacter smithii

# Example 2: Idempotent (applying again changes nothing)
trim_taxonomy_prefix(trimmed)
#>    Kingdom              Genus                    Species
#> 1 Bacteria        Escherichia           Escherichia coli
#> 2  Archaea Methanobrevibacter Methanobrevibacter smithii

# Example 3: Basic usage with phyloseq objects
data("GlobalPatterns", package = "phyloseq")
GlobalPatterns <- add_taxonomy_prefix(GlobalPatterns)
head(phyloseq::tax_table(GlobalPatterns)[, c("Kingdom", "Phylum", "Order")])
#> Taxonomy Table:     [6 taxa by 3 taxonomic ranks]:
#>        Kingdom      Phylum             Order            
#> 549322 "k__Archaea" "p__Crenarchaeota" NA               
#> 522457 "k__Archaea" "p__Crenarchaeota" NA               
#> 951    "k__Archaea" "p__Crenarchaeota" "o__Sulfolobales"
#> 244423 "k__Archaea" "p__Crenarchaeota" NA               
#> 586076 "k__Archaea" "p__Crenarchaeota" NA               
#> 246140 "k__Archaea" "p__Crenarchaeota" NA               

GlobalPatterns <- trim_taxonomy_prefix(GlobalPatterns)
head(phyloseq::tax_table(GlobalPatterns)[, c("Kingdom", "Phylum", "Order")])
#> Taxonomy Table:     [6 taxa by 3 taxonomic ranks]:
#>        Kingdom   Phylum          Order         
#> 549322 "Archaea" "Crenarchaeota" NA            
#> 522457 "Archaea" "Crenarchaeota" NA            
#> 951    "Archaea" "Crenarchaeota" "Sulfolobales"
#> 244423 "Archaea" "Crenarchaeota" NA            
#> 586076 "Archaea" "Crenarchaeota" NA            
#> 246140 "Archaea" "Crenarchaeota" NA            
```
