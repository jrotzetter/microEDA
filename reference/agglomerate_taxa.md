# Agglomerate Taxa by Taxonomic Rank

Collapses a microbiome dataset (`microEDA` or `phyloseq` object) by
aggregating features of the same taxonomy into higher-level taxonomic
groups based on a specified rank. Optionally removes unclassified taxa,
applies total sum scaling (TSS), or adds standard taxonomic rank
prefixes (e.g., `p__`).

## Usage

``` r
agglomerate_taxa(
  me,
  tax_rank,
  rm_missing = FALSE,
  transform = c("None", "TSS"),
  add_prefix = FALSE
)
```

## Arguments

- me:

  A `microEDA` or `phyloseq` object containing abundance and taxonomic
  data.

- tax_rank:

  `Character` string specifying the taxonomic rank for agglomeration
  (e.g., "Genus", "Family"). Must exist in the tax_table of `me`.

- rm_missing:

  `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
  at the specified rank. If `FALSE`, fills missing values by propagating
  the last known ancestor, labeling them as "Unclassified
  *Last_Known_Parent_Clade*" (e.g., "Unclassified Enterobacteriaceae").

- transform:

  `Character.` Transformation to apply to abundances after
  agglomeration. One of `"None"` (no transformation) or `"TSS"` (Total
  Sum Scaling to relative abundance). If `filtered_taxa` are present in
  a `microEDA` object, they will be included in the transformation
  calculation.

- add_prefix:

  `Logical`. If `TRUE`, adds QIIME-style prefixes (e.g., `k__`, `p__`)
  to taxonomic labels.

## Value

Returns an object of the same class as input (`microEDA` or `phyloseq`),
with taxa agglomerated at the specified rank. Preserves sample data,
phylogenetic tree (pruned), and reference sequences if present in `me`.

## Details

This function performs safe agglomeration by grouping taxa based on
their full lineage up to the target rank, reducing the risk of
incorrectly merging distinct clades. A warning is issued if multiple
distinct higher-rank lineages map to the same group at the target level,
indicating potential annotation inconsistencies.

When multiple sequences (ASVs/OTUs) belong to the same taxonomic group
at the target rank, the ID assigned to the agglomerated feature is taken
from the most abundant sequence (summed across samples) within that
group. This should ensure that the dominant biological variant drives
ASV/OTU labeling, enhancing representativeness.

## Note

This implementation is significantly faster than
[`phyloseq::tax_glom`](https://rdrr.io/pkg/phyloseq/man/tax_glom.html),
especially on large datasets, due to efficient use of vectorized
operations. Benchmarks show speedups of 10x to over 100x compared to
`tax_glom`.

## Examples

``` r
# Example with a phyloseq object
data("GlobalPatterns", package = "phyloseq")
agglom <- agglomerate_taxa(GlobalPatterns,
  tax_rank = "Phylum",
  rm_missing = TRUE, transform = "TSS",
  add_prefix = TRUE
)
agglom
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 66 taxa and 26 samples ]
#> sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 66 taxa by 2 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 66 tips and 65 internal nodes ]
```
