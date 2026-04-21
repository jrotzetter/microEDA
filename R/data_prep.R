#' Fill unclassified taxonomic rank using last known parent clade
#'
#' This internal function replaces missing or empty taxonomic values (including `NA`, `""`, `" "`, `"\t"`)
#' with "Unclassified". For the specified target rank, it further relabels "Unclassified" entries
#' as "Unclassified {Last_Known_Parent}" using the last known parent clade.
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#'   (e.g., Kingdom, Phylum, Class) in hierarchical order, and rows represent taxa.
#' @param taxrank A `character` string specifying the target taxonomic rank (e.g., "Genus")
#'   to process and relabel.
#'
#' @return A `matrix` with:
#'   - Missing/empty values replaced by "Unclassified".
#'   - Only columns from the first up to and including `taxrank`.
#'   - "Unclassified" entries in the `taxrank` column relabeled as "Unclassified {Parent}"
#'     where `{Parent}` is the last non-"Unclassified" value in the row.
#'   - Fully unclassified rows remain as "Unclassified".
#'
#' @examples
#' taxa <- data.frame(
#'   Kingdom = c("Kingdom1", "Kingdom2", "Kingdom3", "Kingdom4", "Kingdom5", "Kingdom6", ""),
#'   Phylum = c("Phylum1", "", "Phylum3", "Phylum4", "Phylum5", "", NA),
#'   Class = c("Class1", "Class2", " ", "\t", NA, "Class6", ""),
#'   Order = c("Order1", "Order2", NA, "Order4", NA, "", " ")
#' )
#' rownames(taxa) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5", "OTU6", "OTU7")
#'
#' .fill_unclassified_rank(taxa, "Order")
#'
#' @keywords internal
#' @noRd
.fill_unclassified_rank <- function(taxa, taxrank, silent = TRUE) {
  stopifnot(is.data.frame(taxa) || is.matrix(taxa))
  stopifnot(is.character(taxrank))
  stopifnot(taxrank %in% colnames(taxa))

  # Convert to character matrix for speed
  taxa <- as.matrix(taxa)

  # Replace all NA's or empty strings
  empty_vals <- c(NA, "", " ", "\t")
  taxa[taxa %in% empty_vals] <- "Unclassified"

  # Keep only columns up to and including 'taxrank'
  idx <- which(colnames(taxa) == taxrank)
  taxa <- taxa[, 1:idx, drop = FALSE]

  # Find rows where target rank is "Unclassified"
  is_unclass <- taxa[, idx] == "Unclassified"
  if (!any(is_unclass)) {
    return(taxa)
  }

  # Find last non-"Unclassified" value for each row
  parent_ranks <- taxa[is_unclass, 1:(idx - 1), drop = FALSE]
  last_valid <- vapply(
    seq_len(nrow(parent_ranks)),
    function(i) {
      row <- parent_ranks[i, , drop = FALSE]
      valid <- tail(row[row != "Unclassified"], 1)
      if (length(valid) == 0) "" else valid
    },
    character(1L)
  )

  if (!silent) {
    n_all_unclass <- sum(last_valid == "")
    if (n_all_unclass > 0) {
      warning(n_all_unclass, " rows are fully unclassified up to rank '", taxrank, "'.")
    }
  }

  # Update target column
  taxa[is_unclass, idx] <- ifelse(last_valid == "", "Unclassified",
    paste("Unclassified", last_valid)
  )

  return(taxa)
}


#' Apply Total Sum Scaling (TSS) to numeric columns
#'
#' @param x A `data.frame` or `matrix` containing numeric columns to be scaled.
#' @param presence_threshold `Numeric` value indicating the minimum threshold for
#' presence; values below this are set to 0. Default is 0.
#' @param prop_scale Scaling factor applied after normalization. Default is 100
#' (i.e., convert to percentages).
#'
#' @returns The input object x with numeric columns transformed by TSS.
#' Non-numeric columns remain unchanged.
#' @keywords internal
#' @noRd
.apply_tss <- function(x, presence_threshold = 0, prop_scale = 100) {
  if (any(x < 0, na.rm = TRUE)) {
    warning("Data contains negative entries. Result from applying TSS may not make sense.")
  }

  # Get all numeric columns
  is_numeric <- vapply(x, is.numeric, logical(1L))

  x[is_numeric] <- lapply(x[is_numeric], function(col) {
    col_thresh <- ifelse(col >= presence_threshold, col, 0)

    # Early exit if all values are zero after thresholding to avoid unnecessary
    # computation
    if (all(col_thresh == 0, na.rm = TRUE)) {
      return(rep(0, length(col_thresh)))
    }

    col_norm <- col_thresh / sum(col_thresh, na.rm = TRUE) * prop_scale
    # Replace NaN values resulting from division (e.g., 0/0) with 0
    ifelse(is.nan(col_norm), 0, col_norm)
  })
  return(x)
}
