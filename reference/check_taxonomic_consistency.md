# Check Taxonomic Consistency Across Higher-Level Classifications

Verifies that each taxon at a specified rank has a consistent
higher-level taxonomic path (e.g., all features labeled as the same
genus belong to the same family, order, etc.). Useful for identifying
inconsistencies (from, e.g., different database versions) in taxonomic
annotations within microbiome datasets.

## Usage

``` r
check_taxonomic_consistency(me, tax_rank = NULL, detailed_report = TRUE)
```

## Arguments

- me:

  A `microEDA` or `phyloseq` object, or a `data.frame/matrix` where rows
  represent taxa and columns represent taxonomic ranks (e.g., Kingdom,
  Phylum, Class, etc.).

- tax_rank:

  `Character` string specifying the taxonomic rank at which to check
  consistency (e.g., "Genus", "Family"). If `NULL`, defaults to the
  lowest (last) taxonomic rank present in the data.

- detailed_report:

  `Logical`. If `TRUE`, provides a detailed report showing at which
  higher taxonomic levels inconsistencies occur. If `FALSE`, only
  returns the names of inconsistent taxa.

## Value

Invisibly returns one of the following:

- `NULL` if no inconsistencies are found.

- A `character` vector of inconsistent taxon names if
  `detailed_report = FALSE`.

- A `tibble` with detailed per-rank inconsistency flags if
  `detailed_report = TRUE` and inconsistencies exist.

To use the return value, assign the function call to a variable (e.g.,
`result <- check_taxonomic_consistency(me, detailed_report = TRUE)`).

## Details

The function constructs a "clade lineage" by concatenating all taxonomic
ranks up to the specified `tax_rank`. It then groups by the `tax_rank`
value and checks whether multiple clade paths exist for the same taxon
at that rank. If so, a warning is issued.

When `detailed_report = TRUE`, the function performs additional analysis
to identify exactly which higher ranks differ for each inconsistent
taxon.

## Examples

``` r
# Example using a simple data frame
tax_data <- data.frame(
  Kingdom = rep("Bacteria", 4),
  Phylum = c("Firmicutes", "Firmicutes", "Proteobacteria", "Firmicutes"),
  Class = c("Bacilli", "Bacilli", "Gammaproteobacteria", "Bacilli"),
  Genus = c("Lactobacillus", "Lactobacillus", "Escherichia", "Lactobacillus")
)

# Check consistency at Genus level
check_taxonomic_consistency(tax_data, tax_rank = "Genus")
#> No taxonomic inconsistencies detected at rank 'Genus'.
```
