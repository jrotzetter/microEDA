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
