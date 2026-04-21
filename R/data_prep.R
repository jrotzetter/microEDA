#' Fill unclassified taxonomic ranks with parent-level context (specified rank only)
#'
#' This internal function replaces missing or empty values
#' (including `NA`, `""`, `" "`, `"\t"`) with "Unclassified", while further replacing
#' "Unclassified" entries for the specified rank with the name of the last known
#' taxonomic parent.
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#' (e.g., Kingdom, Phylum) and rows represent taxa.
#' @param taxrank `Character` string specifying the target taxonomic rank
#' (e.g., "Genus") to process and relabel.
#'
#' @return A `data.frame` with:
#'   - Missing/empty values replaced by "Unclassified".
#'   - Only columns from the first up to and including `taxrank`.
#'   - "Unclassified" entries in the `taxrank` column relabeled as "Unclassified {Last_Known_Clade}"
#'     (e.g., "Unclassified Firmicutes" if parent is "Firmicutes").
#'   - Handles cascading unclassified levels by searching upward for the nearest valid taxon.
#'
#' @examples
#' taxa <- data.frame(
#'   Kingdom = c("Kingdom1", "Kingdom2", "Kingdom3", "Kingdom4", "Kingdom5"),
#'   Phylum = c("Phylum1", "", "Phylum3", "Phylum4", "Phylum5"),
#'   Class = c("Class1", "Class2", " ", "\t", NA),
#'   Order = c("Order1", "Order2", NA, "Order4", NA)
#' )
#' rownames(taxa) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5")
#'
#' .fill_unclassified_rank(taxa, "Order")
#'
#' @keywords internal
#' @noRd
.fill_unclassified_rank <- function(taxa, taxrank) {
  stopifnot(is.matrix(taxa) || is.data.frame(taxa))
  stopifnot(is.character(taxrank))

  # Replace all NA's or empty strings
  empty_vals <- c(NA, "", " ", "\t")

  taxa[] <- apply(taxa, 2, function(col) {
    col <- as.character(col)
    col[is.na(col) | col %in% c("", " ", "\t")] <- "Unclassified"
    col
  })

  # Keep only columns up to and including 'taxrank'
  idx <- which(colnames(taxa) == taxrank)
  taxa <- taxa[, 1:idx, drop = FALSE]

  # Replace "Unclassified" with "Unclassified [previous_level]"
  is_unclass <- taxa[, taxrank] == "Unclassified"
  if (any(is_unclass)) {
    prev_col_idx <- which(colnames(taxa) == taxrank) - 1
    taxa[, taxrank][is_unclass] <- paste0(
      "Unclassified ",
      as.character(taxa[, prev_col_idx][is_unclass])
    )
  }

  # Handle multiple "Unclassified" levels when previous level was also unclassified
  for (i in 1:(ncol(taxa) - 1)) {
    is_double <- taxa[, taxrank] == "Unclassified Unclassified"
    if (any(is_double)) {
      col_to_fetch <- max(1, prev_col_idx - i + 1) # Avoid negative index
      taxa[, taxrank][is_double] <- paste0(
        "Unclassified ",
        as.character(taxa[, col_to_fetch][is_double])
      )
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
