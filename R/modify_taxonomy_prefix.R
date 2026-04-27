#' Add taxonomic prefixes to tax_table columns
#'
#' Internal function to add QIIME-style two-underscore prefixes (e.g., "k__", "p__")
#' to values in a taxonomy matrix if they are not already present. The function
#' is idempotent, meaning it will not add duplicate prefixes.
#'
#' @param tax_mat A `character` matrix or data.frame where column names correspond
#'  to taxonomic ranks (e.g., "Kingdom", "Phylum"). Values are taxonomic labels
#'  to be prefixed.
#' @param prefix_map A named `character` vector mapping taxonomic ranks
#'  (e.g., "Kingdom") to their prefixes (e.g., "k__"). Defaults to QIIME-style prefixes.
#' @return Returns the input matrix/data.frame with prefixes added to appropriate
#'  cells in relevant columns. Only columns matching taxonomic ranks in `prefix_map`
#'  and values not already prefixed are modified. Empty, blank, or `NA` values
#'  are left unmodified.
#' @details
#'  - Applies prefixes only to columns present in both the input `tax_mat` and
#'  the prefix_map ("Kingdom", "Phylum", etc.).
#'  - Uses a single regex pattern to check if values already start with any
#'  defined prefix, preventing redundant prefixing.
#'
#' @seealso [trim_taxonomy_prefix()] for removing prefixes.
#'
#' @examples
#' tax_data <- data.frame(
#'   Kingdom = c("Bacteria", "Archaea"),
#'   Genus = c("Escherichia", "Methanobrevibacter"),
#'   Species = c("Escherichia coli", "Methanobrevibacter smithii"),
#'   stringsAsFactors = FALSE
#' )
#' prefixed <- add_taxonomy_prefix(tax_data)
#' print(prefixed)
#'
#' # Idempotent: applying again changes nothing
#' add_taxonomy_prefix(prefixed)
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stringi stri_detect_regex
.add_taxonomy_prefix <- function(tax_mat, prefix_map = c(
                                   "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                                   "Order" = "o__", "Family" = "f__", "Genus" = "g__",
                                   "Species" = "s__", "Strain" = "t__"
                                 )) {
  # Input validation
  if (!is.matrix(tax_mat) && !is.data.frame(tax_mat)) {
    stop("`tax_mat` must be a character matrix or data.frame.")
  }
  if (is.matrix(tax_mat) && !is.character(tax_mat)) {
    stop("`tax_mat` must be of character type.")
  }

  # Find which ranks are actually present in the matrix
  common_ranks <- intersect(names(prefix_map), colnames(tax_mat))

  if (length(common_ranks) == 0) {
    stop("No matching taxonomic ranks found between `prefix_map` and `tax_mat`.")
  }

  # Create a dynamic regex pattern from all prefixes
  all_prefixes_pattern <- paste0("^(", paste(prefix_map, collapse = "|"), ")")

  # Apply prefix only to present ranks
  for (rank in common_ranks) {
    # Check for values that are empty or whitespace-only
    is_blank <- !nzchar(trimws(tax_mat[, rank]))

    # Check if the value does NOT already start with any of the defined prefixes
    # needs_prefix <- !grepl(all_prefixes_pattern, tax_mat[, rank], perl = TRUE)

    needs_prefix <- !stringi::stri_detect_regex(tax_mat[, rank], all_prefixes_pattern)
    # stringi::stri_detect_regex() returns NA for NA inputs, which propagates
    # through the logical indexing. Need to explicitly convert NA values in the
    # logical vector to FALSE, ensuring they won't trigger prefix addition, plus
    # ensure that empty or blank entries are not prefixed
    needs_prefix[is.na(needs_prefix) | is_blank] <- FALSE

    if (any(needs_prefix)) {
      # Add the rank-specific prefix
      tax_mat[needs_prefix, rank] <- paste0(prefix_map[rank], tax_mat[needs_prefix, rank])
    }
  }
  return(tax_mat)
}


#' Add taxonomic prefixes to tax_table columns
#'
#' Adds QIIME-style two-underscore prefixes (e.g., "k__", "p__") to values in a
#' taxonomy matrix or data.frame if they are not already present. The function
#' is idempotent, meaning it will not add duplicate prefixes when reapplied.
#'
#' This function supports both character matrices/data.frames and `microEDA` or
#' `phyloseq` objects.
#' It only modifies columns matching known taxonomic ranks (as defined in `prefix_map`),
#' leaving other columns and missing/empty values (`NA`, "", " ") unchanged.
#'
#' @param me A `character` matrix, `data.frame`, `microEDA` or `phyloseq`-class
#'   object containing a `tax_table`. If a `microEDA` or `phyloseq` object is passed,
#'   its internal taxonomyTable is modified.
#' @param prefix_map A named `character` vector mapping taxonomic rank names
#'   (e.g., "Kingdom", "Phylum") to their corresponding prefixes (e.g., "k__", "p__").
#'   Defaults to standard QIIME-style prefixes for common ranks including "Strain" ("t__").
#'
#' @return Returns the input object with prefixes added to appropriate cells:
#'   - For matrices/data.frames: returns the modified matrix or data.frame.
#'   - For `microEDA`/`phyloseq` objects: returns the updated `microEDA`/`phyloseq`
#'    object with prefixed `tax_table`.
#'
#' @details
#'   - Only columns in `me` that match keys in `prefix_map` are processed.
#'   - Values already starting with any valid prefix (e.g., "k__Bacteria") are skipped.
#'   - Empty strings (`""`), `NA`, and whitespace-only entries are preserved without modification.
#'
#' @seealso
#'   [trim_taxonomy_prefix()] to remove prefixes.
#'
#'   [tax_table()] to access taxonomy data.
#'
#' @examples
#' # Example 1: Basic usage with data.frame
#' tax_data <- data.frame(
#'   Kingdom = c("Bacteria", "Archaea"),
#'   Phylum = c("Proteobacteria", "Euryarchaeota"),
#'   Genus = c("Escherichia", "Methanobrevibacter"),
#'   Species = c("Escherichia coli", "Methanobrevibacter smithii"),
#'   stringsAsFactors = FALSE
#' )
#'
#' prefixed <- add_taxonomy_prefix(tax_data)
#' prefixed
#'
#' # Example 2: Idempotent — applying again changes nothing
#' add_taxonomy_prefix(prefixed)
#'
#' # Example 3: Basic usage with phyloseq objects
#' data("GlobalPatterns", package = "phyloseq")
#' GlobalPatterns <- add_taxonomy_prefix(GlobalPatterns)
#' head(phyloseq::tax_table(GlobalPatterns)[, c("Kingdom", "Phylum", "Order")])
#'
#' @export
add_taxonomy_prefix <- function(me, prefix_map = c(
                                  "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                                  "Order" = "o__", "Family" = "f__", "Genus" = "g__",
                                  "Species" = "s__", "Strain" = "t__"
                                )) {
  # Input validation
  if (inherits(me, "microEDA") || inherits(me, "phyloseq")) {
    # Proceed with microEDA or phyloseq processing
    # slot(me, "tax_table") <- .add_taxonomy_prefix(slot(me, "tax_table"), prefix_map = prefix_map)
    tax_table(me) <- .add_taxonomy_prefix(tax_table(me), prefix_map = prefix_map)
  } else if (is.matrix(me) || is.data.frame(me)) {
    if (is.matrix(me) && !is.character(me)) {
      stop("`me` must be of character type if it's a matrix.")
    }
    # Proceed with matrix/data.frame processing
    me <- .add_taxonomy_prefix(me, prefix_map = prefix_map)
  } else {
    stop("'me' must be a microEDA, phyloseq, or a character matrix/data.frame.")
  }
  return(me)
}


#' Remove taxonomic prefixes from tax_table columns
#'
#' Internal function to remove QIIME-style two-underscore prefixes (e.g., "k__", "p__")
#' from values in a taxonomy matrix if present. The function is idempotent,
#' meaning repeated application has no additional effect.
#'
#' @param tax_mat A `character` matrix or data.frame where column names correspond
#'  to taxonomic ranks (e.g., "Kingdom", "Phylum"). Values may contain prefixed labels.
#' @param prefix_map A named `character` vector mapping taxonomic ranks
#'  (e.g., "Kingdom") to their prefixes (e.g., "k__"). Defaults to QIIME-style prefixes.
#' @return Returns the input matrix/data.frame with prefixes removed from appropriate
#'  cells in relevant columns. Only values in matching columns and starting with the
#'  correct prefix are modified. Empty, blank, or `NA` values are left unmodified.
#' @details
#'  - Applies only to columns present in both `tax_mat` and `prefix_map`.
#'  - Uses exact prefix-to-rank mapping: e.g., only `k__` is stripped from "Kingdom".
#'  - Relies on `stringi::stri_detect_regex` for robust pattern matching.
#'
#' @seealso [add_taxonomy_prefix()] for adding prefixes.
#'
#' @examples
#' tax_data_with_prefix <- data.frame(
#'   Kingdom = c("k__Bacteria", "k__Archaea"),
#'   Genus = c("g__Escherichia", "g__Methanobrevibacter"),
#'   Species = c("s__Escherichia coli", "s__Methanobrevibacter smithii"),
#'   stringsAsFactors = FALSE
#' )
#' trimmed <- trim_taxonomy_prefix(tax_data_with_prefix)
#' print(trimmed)
#'
#' # Idempotent: applying again changes nothing
#' trim_taxonomy_prefix(trimmed)
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stringi stri_detect_regex
#' @importFrom stringr str_escape
.trim_taxonomy_prefix <- function(tax_mat, prefix_map = c(
                                    "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                                    "Order" = "o__", "Family" = "f__", "Genus" = "g__",
                                    "Species" = "s__", "Strain" = "t__"
                                  )) {
  # Input validation
  if (!is.matrix(tax_mat) && !is.data.frame(tax_mat)) {
    stop("`tax_mat` must be a character matrix or data.frame.")
  }
  if (is.matrix(tax_mat) && !is.character(tax_mat)) {
    stop("`tax_mat` must be of character type.")
  }

  # Find which ranks are actually present in the matrix
  common_ranks <- intersect(names(prefix_map), colnames(tax_mat))

  if (length(common_ranks) == 0) {
    stop("No matching taxonomic ranks found between `prefix_map` and `tax_mat`.")
  }

  # Apply to each relevant column
  for (rank in common_ranks) {
    prefix <- prefix_map[rank]
    # Escape special regex characters in prefix (though "__" is safe, need to
    # consider user provided strings through prefix_map)
    escaped_prefix <- paste0("^", stringr::str_escape(prefix))

    # Identify non-blank, non-NA entries
    is_blank <- !nzchar(trimws(tax_mat[, rank]))
    has_prefix <- stringi::stri_detect_regex(tax_mat[, rank], escaped_prefix)
    # Combine conditions: must have prefix and not be blank/NA
    can_strip <- has_prefix & !is_blank
    # Handle NA in logical vector
    can_strip[is.na(can_strip)] <- FALSE

    if (any(can_strip)) {
      # Remove exactly the prefix at the start
      tax_mat[can_strip, rank] <- substring(tax_mat[can_strip, rank], nchar(prefix) + 1)
    }
  }
  return(tax_mat)
}


#' Remove taxonomic prefixes from tax_table columns
#'
#' Removes QIIME-style two-underscore prefixes (e.g., "k__", "p__") from values
#' in a taxonomy matrix or data.frame if present. The function is idempotent:
#' repeated application has no additional effect.
#'
#' Supports character matrices, data.frames, `microEDA` and `phyloseq` objects.
#' Only modifies cells in columns matching the keys of `prefix_map` and only when
#' the value starts with the corresponding prefix. Missing (`NA`), empty, or
#' unprefixed values are left unchanged.
#'
#' @param me A `character` matrix, `data.frame`, `microEDA` or a `phyloseq`-class
#'   object containing a `tax_table`. If a `microEDA` or `phyloseq` object is
#'   passed, its taxonomyTable is modified.
#' @param prefix_map A named `character` vector mapping taxonomic rank names
#'   (e.g., "Kingdom", "Phylum") to their corresponding prefixes (e.g., "k__", "p__").
#'   Defaults to standard QIIME-style prefixes including "Strain" ("t__").
#'
#' @return Returns the input object with prefixes removed:
#'   - For matrices/data.frames: returns a modified `character` matrix or data.frame.
#'   - For `microEDA`/`phyloseq` objects: returns the updated `microEDA`/`phyloseq` object with cleaned `tax_table`.
#'
#' @details
#'   - Only columns in `me` that match names in `prefix_map` are processed.
#'   - Prefix removal is exact: e.g., "k__" is only removed from "Kingdom" column if it starts the value.
#'   - Uses efficient string processing (e.g., via `stringi`) to detect and remove prefixes.
#'   - Preserves `NA`, `""`, and whitespace-only entries.
#'
#' @seealso
#'   [add_taxonomy_prefix()] to add prefixes.
#'
#'   [tax_table()] to access taxonomy data.
#'
#' @examples
#' # Example 1: Basic usage with data.frame
#' prefixed_tax_data <- data.frame(
#'   Kingdom = c("k__Bacteria", "k__Archaea"),
#'   Genus = c("g__Escherichia", "g__Methanobrevibacter"),
#'   Species = c("s__Escherichia coli", "s__Methanobrevibacter smithii"),
#'   stringsAsFactors = FALSE
#' )
#' prefixed_tax_data
#' trimmed <- trim_taxonomy_prefix(prefixed_tax_data)
#' trimmed
#'
#' # Example 2: Idempotent — applying again changes nothing
#' trim_taxonomy_prefix(trimmed)
#'
#' # Example 3: Basic usage with phyloseq objects
#' data("GlobalPatterns", package = "phyloseq")
#' GlobalPatterns <- add_taxonomy_prefix(GlobalPatterns)
#' head(phyloseq::tax_table(GlobalPatterns)[, c("Kingdom", "Phylum", "Order")])
#'
#' GlobalPatterns <- trim_taxonomy_prefix(GlobalPatterns)
#' head(phyloseq::tax_table(GlobalPatterns)[, c("Kingdom", "Phylum", "Order")])
#'
#' @export
trim_taxonomy_prefix <- function(me, prefix_map = c(
                                   "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
                                   "Order" = "o__", "Family" = "f__", "Genus" = "g__",
                                   "Species" = "s__", "Strain" = "t__"
                                 )) {
  # Input validation
  if (inherits(me, "microEDA") || inherits(me, "phyloseq")) {
    # Proceed with microEDA or phyloseq processing
    tax_table(me) <- .trim_taxonomy_prefix(tax_table(me), prefix_map = prefix_map)
  } else if (is.matrix(me) || is.data.frame(me)) {
    if (is.matrix(me) && !is.character(me)) {
      stop("`me` must be of character type if it's a matrix.")
    }
    # Proceed with matrix/data.frame processing
    me <- .trim_taxonomy_prefix(me, prefix_map = prefix_map)
  } else {
    stop("'me' must be a microEDA, phyloseq, or a character matrix/data.frame.")
  }
  return(me)
}
