# Add taxonomic prefixes to tax_table columns

Adds QIIME-style two-underscore prefixes (e.g., "k\_\_", "p\_\_") to
values in a taxonomy matrix or data.frame if they are not already
present. The function is idempotent, meaning it will not add duplicate
prefixes when reapplied.

## Usage

``` r
add_taxonomy_prefix(
  me,
  prefix_map = c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__", Family =
    "f__", Genus = "g__", Species = "s__", Strain = "t__")
)
```

## Arguments

- me:

  A `character` matrix, `data.frame`, `microEDA` or `phyloseq`-class
  object containing a `tax_table`. If a `microEDA` or `phyloseq` object
  is passed, its internal taxonomyTable is modified.

- prefix_map:

  A named `character` vector mapping taxonomic rank names (e.g.,
  "Kingdom", "Phylum") to their corresponding prefixes (e.g., "k\_\_",
  "p\_\_"). Defaults to standard QIIME-style prefixes for common ranks
  including "Strain" ("t\_\_").

## Value

Returns the input object with prefixes added to appropriate cells:

- For matrices/data.frames: returns the modified matrix or data.frame.

- For `microEDA`/`phyloseq` objects: returns the updated
  `microEDA`/`phyloseq` object with prefixed `tax_table`.

## Details

This function supports both character matrices/data.frames and
`microEDA` or `phyloseq` objects. It only modifies columns matching
known taxonomic ranks (as defined in `prefix_map`), leaving other
columns and missing/empty values (`NA`, "", " ") unchanged.

- Only columns in `me` that match keys in `prefix_map` are processed.

- Values already starting with any valid prefix (e.g., "k\_\_Bacteria")
  are skipped.

- Empty strings (`""`), `NA`, and whitespace-only entries are preserved
  without modification.

## See also

[`trim_taxonomy_prefix()`](https://jrotzetter.github.io/microEDA/reference/trim_taxonomy_prefix.md)
to remove prefixes.

[`tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)
to access taxonomy data.

## Examples

``` r
# Example 1: Basic usage with data.frame
tax_data <- data.frame(
  Kingdom = c("Bacteria", "Archaea"),
  Phylum = c("Proteobacteria", "Euryarchaeota"),
  Genus = c("Escherichia", "Methanobrevibacter"),
  Species = c("Escherichia coli", "Methanobrevibacter smithii"),
  stringsAsFactors = FALSE
)

prefixed <- add_taxonomy_prefix(tax_data)
prefixed
#>       Kingdom            Phylum                 Genus
#> 1 k__Bacteria p__Proteobacteria        g__Escherichia
#> 2  k__Archaea  p__Euryarchaeota g__Methanobrevibacter
#>                         Species
#> 1           s__Escherichia coli
#> 2 s__Methanobrevibacter smithii

# Example 2: Idempotent (applying again changes nothing)
add_taxonomy_prefix(prefixed)
#>       Kingdom            Phylum                 Genus
#> 1 k__Bacteria p__Proteobacteria        g__Escherichia
#> 2  k__Archaea  p__Euryarchaeota g__Methanobrevibacter
#>                         Species
#> 1           s__Escherichia coli
#> 2 s__Methanobrevibacter smithii

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
```
