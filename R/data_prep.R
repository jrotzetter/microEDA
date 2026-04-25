#' Fill unclassified taxonomic rank using last known parent clade
#'
#' This internal function replaces missing or empty taxonomic entries (including `NA`, `""`, `" "`, `"\t"`)
#' with "Unclassified". For the specified target rank, it further relabels "Unclassified" entries
#' as "Unclassified {Last_Known_Parent}" using the last known parent clade.
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#'   (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species) in hierarchical order,
#'   and rows represent taxa.
#' @param taxrank A `character` string specifying the target taxonomic rank (e.g., "Genus")
#'   to process and relabel.
#' @param silent A `logical` value indicating whether to suppress warnings.
#'   If `FALSE`, a warning is issued if any rows are fully unclassified across all ranks.
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
      warning(n_all_unclass, " row(s) are fully unclassified up to rank '", taxrank, "'.")
    }
  }

  # Update target column
  taxa[is_unclass, idx] <- ifelse(last_valid == "", "Unclassified",
    paste("Unclassified", last_valid)
  )

  return(taxa)
}


#' Fill all unclassified taxonomic ranks using the last known parent clade
#'
#' This internal function replaces missing or empty taxonomic entries (including `NA`, `""`, `" "`, `"\t"`)
#' with "Unclassified". It then processes all taxonomic ranks sequentially, relabeling "Unclassified"
#' entries as "Unclassified {Last_Known_Parent}" using the last known valid parent clade from higher ranks.
#'
#' The function propagates the most recent non-unclassified taxon name across lower ranks,
#' ensuring consistent hierarchical context for unclassified taxa. Fully unclassified rows
#' (i.e., all ranks are "Unclassified") remain labeled as "Unclassified".
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#'   (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species) in hierarchical order,
#'   and rows represent taxa.
#' @param silent A `logical` value indicating whether to suppress warnings.
#'   If `FALSE`, a warning is issued if any rows are fully unclassified across all ranks.
#'
#' @return A `matrix` of the same dimensions as input, with:
#'   - All missing/empty values replaced by "Unclassified".
#'   - "Unclassified" entries in each rank relabeled as "Unclassified {Parent}"
#'     where `{Parent}` is the last non-unclassified taxon in the row up to that point.
#'   - Fully unclassified rows remain as "Unclassified" in all columns.
#'   - All original columns preserved (no column subsetting).
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
#' .fill_unclassified_all(taxa, silent = FALSE)
#'
#' @keywords internal
#' @noRd
.fill_unclassified_all <- function(taxa, silent = TRUE) {
  stopifnot(is.data.frame(taxa) || is.matrix(taxa))

  # Convert to character matrix for speed
  taxa <- as.matrix(taxa)

  # Replace all NA's or empty strings
  empty_vals <- c(NA, "", " ", "\t")
  taxa[taxa %in% empty_vals] <- "Unclassified"

  # Track last valid (non-unclassified) name per row across ranks
  last_valid <- rep("", nrow(taxa))

  for (j in seq_len(ncol(taxa))) {
    col <- taxa[, j]
    is_unclass <- col == "Unclassified"
    # Replace unclassified with "Unclassified <Parent>" using last valid name
    col[is_unclass] <- ifelse(last_valid[is_unclass] == "",
      "Unclassified",
      paste("Unclassified", last_valid[is_unclass])
    )
    taxa[, j] <- col

    # Update last_valid with names of "true" taxa (not "Unclassified" or
    # "Unclassified <Parent>") so they propagate to lower ranks
    true_name <- !grepl("^Unclassified($| )", taxa[, j])
    last_valid[true_name] <- taxa[true_name, j]
  }

  if (!silent) {
    all_unclass <- apply(taxa, 1, function(x) all(grepl("^Unclassified($| )", x)))
    n_all_unclass <- sum(all_unclass)
    if (n_all_unclass > 0) {
      warning(n_all_unclass, " row(s) are fully unclassified across all ranks.")
    }
  }

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

  # Handle data.frame and matrix differently
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("Matrix must be numeric.")
    # # Vectorized thresholding
    # x[x < presence_threshold] <- 0
    # col_sums <- colSums(x, na.rm = TRUE)
    # # Replace zero sums to avoid NaN
    # col_sums[col_sums == 0] <- 1
    # x <- sweep(x, 2, col_sums, "/") * prop_scale
    # x[is.nan(x)] <- 0

    # For speed, convert to data.frame, process, then convert back as data.frames
    # store columns as list elements, allowing fast column access (matrices are
    # vectors with dimensions; column extraction requires more computation)
    x_df <- as.data.frame(x)
    x_df[] <- lapply(x_df, function(col) {
      col_thresh <- ifelse(col >= presence_threshold, col, 0)
      if (all(col_thresh == 0, na.rm = TRUE)) {
        return(rep(0, length(col)))
      }
      col_norm <- col_thresh / sum(col_thresh, na.rm = TRUE) * prop_scale
      ifelse(is.nan(col_norm), 0, col_norm)
    })
    x <- as.matrix(x_df)
  } else {
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
  }
  return(x)
}
