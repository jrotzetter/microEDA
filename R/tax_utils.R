#' Taxonomic rank abbreviations
#'
#' Named character vector mapping full taxonomic rank names to their abbreviations.
#' Used internally by `.trim_rank_prefix()` and `.add_rank_prefix()`.
#' @keywords internal
#' @noRd
.tax_ranks <- c(
  Kingdom = "k",
  Phylum = "p",
  Class = "c",
  Order = "o",
  Family = "f",
  Genus = "g",
  Species = "s",
  Strain = "t"
)


#' Remove taxonomic rank prefix
#'
#' Internal function to remove the rank prefix (e.g. "g__") from a taxonomic
#' classification string.
#'
#' @param x A `character` vector containing taxonomic strings with rank prefixes.
#' @return A `character` vector with rank prefixes removed.
#' @keywords internal
#' @noRd
.trim_rank_prefix <- function(x) {
  pattern <- paste(.tax_ranks, collapse = "")
  pattern <- paste0("^[", pattern, "]__")

  match <- gsub(pattern, "", as.character(x))
  return(match)
}


#' Add taxonomic rank prefixes
#'
#' Internal function to add rank prefixes (e.g. "g__") to columns in a
#' taxonomic table.
#'
#' @param tax_table A data frame or tibble containing taxonomic classifications.
#' @param cols Character vector specifying the column names in `tax_table` to
#'   which prefixes should be added.
#' @return The input `tax_table` with rank prefixes added to the specified columns.
#' @keywords internal
#' @noRd
.add_rank_prefix <- function(tax_table, cols) {
  for (tax_rank in cols) {
    prefix <- .tax_ranks[tax_rank]
    tax_table[[tax_rank]] <- paste0(prefix, "__", tax_table[[tax_rank]])
  }
  return(tax_table)
}


#' Map abbreviated or lowercase taxonomic rank to full name
#'
#' Converts short-form taxonomic rank codes (e.g., "p", "g") or lowercase names
#' to their full names (e.g., "Phylum", "Genus").
#'
#' @param tax_rank `Character` vector with abbreviated or lowercase taxonomic rank name.
#' @return `Character` vector with the full rank name (e.g., "Family").
#' @keywords internal
#' @noRd
.get_full_tax_rank <- function(tax_rank) {
  lookup <- c(
    "k" = "Kingdom",
    "p" = "Phylum",
    "c" = "Class",
    "o" = "Order",
    "f" = "Family",
    "g" = "Genus",
    "s" = "Species",
    "t" = "Strain",
    "kingdom" = "Kingdom",
    "phylum" = "Phylum",
    "class" = "Class",
    "order" = "Order",
    "family" = "Family",
    "genus" = "Genus",
    "species" = "Species",
    "strain" = "Strain"
  )
  tax_rank <- unname(lookup[tolower(tax_rank)])
}
