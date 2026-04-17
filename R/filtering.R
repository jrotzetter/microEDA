#' Filter Features Based on Abundance and Prevalence
#'
#' An internal function to filter features (i.e., taxa) from a microbiome dataset.
#'
#' @param x A `data.frame` containing feature abundance data, with samples as
#' columns and taxa names in a `tax_ID` column.
#' @param metadata A `data.frame` containing sample metadata; required if
#'  `group_var` is specified.
#' @param group_var A `character` string specifying a column in `metadata` to
#' group samples by for stratified filtering.
#' @param min_abundance A non-negative `numeric` value; features must have
#' abundance >= to this value to pass.
#' @param min_prevalence A non-negative `numeric` value; interpreted as an
#' absolute sample count if >=1, or as a proportion of samples if <1.
#' @param abundance_criterion A `character` string; either `"prevalence"`
#' (default) to filter based on how often a feature meets `min_abundance`, or
#' `"mean"` to also require the mean abundance passing the threshold.
#' @param group_requirement A `character` string; `"any"` (default) if a feature
#' must pass in at least one group, or `"all"` if it must pass in all groups.
#' @param keep_other A `logical` value; if `TRUE` and filtering is applied, adds
#' an "Other" row summarizing the abundance of filtered features.
#'
#' @return A filtered `data.frame` with features that meet the criteria,
#' optionally including an "Other" row. A message reports the number and names
#' of filtered features.
#'
#' @examples
#' df <- data.frame(
#'   tax_ID = c("tax1", "tax2", "tax3", "tax4"),
#'   SampleA = c(10, 30, 60, 0),
#'   SampleB = c(20, 0, 80, 0),
#'   SampleC = c(0, 40, 60, 0),
#'   SampleD = c(0, 0, 0, 0)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("SampleA", "SampleB", "SampleC", "SampleD"),
#'   Category = c("Cat1", "Cat1", "Cat2", "Cat3")
#' )
#'
#' .filter_features(df,
#'   metadata = metadata, group_var = "Category",
#'   min_abundance = 10, min_prevalence = 2
#' )
#'
#' @keywords internal
#' @noRd
.filter_features <- function(x,
                             metadata = NULL,
                             group_var = NULL,
                             min_abundance = 0,
                             min_prevalence = 0,
                             abundance_criterion = c("prevalence", "mean"),
                             group_requirement = c("any", "all"),
                             keep_other = TRUE) {
  abundance_criterion <- match.arg(abundance_criterion)
  group_requirement <- match.arg(group_requirement)

  # Validate
  if (!is.numeric(min_abundance) || min_abundance < 0) {
    stop("min_abundance must be non-negative numeric")
  }

  if (!is.numeric(min_prevalence) || min_prevalence < 0) {
    stop("min_prevalence must be non-negative numeric")
  }

  is_numeric <- vapply(x, is.numeric, logical(1))
  abund_data <- x[is_numeric]
  tax_ids <- x$tax_ID

  # No grouping
  if (is.null(group_var)) {
    total_samples <- ncol(abund_data)

    # Interpret min_prevalence: if >=1, treat as sample count, else as proportion
    min_samples <- if (min_prevalence >= 1) {
      as.integer(min_prevalence)
    } else {
      ceiling(total_samples * min_prevalence)
    }

    if (abundance_criterion == "prevalence") {
      above_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(above_threshold, na.rm = TRUE)
      keep <- prevalence >= min_samples
    } else { # "mean"
      avg_abund <- rowMeans(abund_data, na.rm = TRUE)
      above_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(above_threshold, na.rm = TRUE)
      keep <- avg_abund >= min_abundance & prevalence >= min_samples
    }
    filtered_taxa <- tax_ids[keep]
  } else {
    # Grouped filtering
    df_long <- tidyr::pivot_longer(x,
      cols = all_of(names(abund_data)),
      names_to = "Sample", values_to = "Abundance"
    ) |>
      dplyr::inner_join(metadata, by = "Sample")

    filter_fn <- if (group_requirement == "all") all else any

    if (abundance_criterion == "prevalence") {
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], tax_ID) |>
        dplyr::summarise(
          pass = sum(Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(n_samples * min_prevalence))) |>
        dplyr::group_by(tax_ID) |>
        dplyr::filter(filter_fn(pass >= min_n)) |>
        dplyr::pull(tax_ID)
    } else { # "mean"
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], tax_ID) |>
        dplyr::summarise(
          avg_abund = mean(Abundance, na.rm = TRUE),
          pass = sum(Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(n_samples * min_prevalence))) |>
        dplyr::group_by(tax_ID) |>
        dplyr::filter(filter_fn(avg_abund >= min_abundance & pass >= min_n)) |>
        dplyr::pull(tax_ID)
    }
    filtered_taxa <- unique(keep)
  }

  filtered_data <- x[x$tax_ID %in% filtered_taxa, , drop = FALSE]

  # Add "Other" row if taxa were likely removed
  if (keep_other && (min_abundance > 0 || min_prevalence > 0)) {
    initial_colsums <- colSums(x[is_numeric])
    filtered_data <- .add_other_row(filtered_data, initial_colsums, tax_col = "tax_ID")
  }

  # Report on removed taxa
  n_filtered <- nrow(x) - nrow(filtered_data[filtered_data$tax_ID != "Other", ])
  if (n_filtered > 0) {
    message("Total filtered features: ", n_filtered, "\n")
    dropped <- setdiff(x$tax_ID, filtered_data$tax_ID)
    message("Filtered features: ", toString(dropped), "\n")
  }

  return(filtered_data)
}


#' @keywords internal
#' @noRd
.add_other_row <- function(filtered_df, initial_colsums, tax_col = "tax_ID") {
  # TODO: prevent function from adding an additional 'Other' row if already
  # present and merge instead
  is_numeric <- vapply(filtered_df, is.numeric, logical(1L))

  filtered_colsums <- colSums(filtered_df[, is_numeric], na.rm = TRUE)

  other_row <- as.data.frame(t(initial_colsums - filtered_colsums))
  other_row[[tax_col]] <- "Other"

  df_combined <- dplyr::bind_rows(filtered_df, other_row)

  return(df_combined)
}
