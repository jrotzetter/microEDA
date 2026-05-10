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
#' @param keep_filtered A `logical` value; if `TRUE` and filtering is applied, adds
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
                             keep_filtered = TRUE) {
  abundance_criterion <- match.arg(abundance_criterion)
  group_requirement <- match.arg(group_requirement)

  # Validate
  if (!is.numeric(min_abundance) || min_abundance < 0) {
    stop("min_abundance must be non-negative numeric")
  }

  if (!is.numeric(min_prevalence) || min_prevalence < 0) {
    stop("min_prevalence must be non-negative numeric")
  }

  is_numeric <- vapply(x, is.numeric, logical(1L))
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
      pass_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(pass_threshold, na.rm = TRUE)
      keep <- prevalence >= min_samples
    } else { # "mean"
      avg_abund <- rowMeans(abund_data, na.rm = TRUE)
      pass_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(pass_threshold, na.rm = TRUE)
      keep <- avg_abund >= min_abundance & prevalence >= min_samples
    }
    filtered_taxa <- tax_ids[keep]
  } else {
    # Grouped filtering
    df_long <- tidyr::pivot_longer(x,
      cols = dplyr::all_of(names(abund_data)),
      names_to = "Sample", values_to = "Abundance"
    ) |>
      dplyr::inner_join(metadata, by = "Sample")

    filter_fn <- if (group_requirement == "all") all else any

    if (abundance_criterion == "prevalence") {
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], .data$tax_ID) |>
        dplyr::summarise(
          pass = sum(.data$Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(.data$n_samples * min_prevalence))) |>
        dplyr::group_by(.data$tax_ID) |>
        dplyr::filter(filter_fn(.data$pass >= .data$min_n)) |>
        dplyr::pull(.data$tax_ID)
    } else { # "mean"
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], .data$tax_ID) |>
        dplyr::summarise(
          avg_abund = mean(.data$Abundance, na.rm = TRUE),
          pass = sum(.data$Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(.data$n_samples * min_prevalence))) |>
        dplyr::group_by(.data$tax_ID) |>
        dplyr::filter(filter_fn(.data$avg_abund >= min_abundance & .data$pass >= .data$min_n)) |>
        dplyr::pull(.data$tax_ID)
    }
    filtered_taxa <- unique(keep)
  }

  filtered_data <- x[x$tax_ID %in% filtered_taxa, , drop = FALSE]

  # Add "Other" row if taxa were likely removed
  if (keep_filtered && (min_abundance > 0 || min_prevalence > 0)) {
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
  is_numeric <- vapply(filtered_df, is.numeric, logical(1L))

  filtered_colsums <- colSums(filtered_df[, is_numeric], na.rm = TRUE)

  # Check for existing 'Other' rows
  other_mask <- filtered_df[[tax_col]] == "Other"

  # Sum all 'Other' rows if any exist
  if (any(other_mask)) {
    other_sum <- colSums(filtered_df[other_mask, is_numeric, drop = FALSE], na.rm = TRUE)
    # Remove all existing 'Other' rows
    filtered_df <- filtered_df[!other_mask, ]
  } else {
    other_sum <- 0
  }

  # Compute difference between initial data and filtered data
  # - other_sum: Adjust by removing prior "Other" as we don't want to exclude
  # them from other_vals
  other_vals <- initial_colsums - (filtered_colsums - other_sum)

  # Convert vector to one row data.frame and assign tax_ID
  other_row <- as.data.frame(t(other_vals))
  other_row[[tax_col]] <- "Other"


  df_combined <- dplyr::bind_rows(filtered_df, other_row)

  return(df_combined)
}


#' Filter Features Based on Abundance and Prevalence
#'
#' Filters taxa (features) from a `microEDA` or `phyloseq` object based on minimum abundance
#' and prevalence thresholds. Filtering can be applied globally or within groups
#' defined by a metadata variable. Optionally adds an "Other" category summarizing
#' filtered features.
#'
#' @param me A `microEDA` or `phyloseq` object containing microbiome abundance data and associated metadata.
#' @param group_var (Optional) Name of the sample metadata variable to define groups
#'   for stratified filtering. If `NULL`, filtering is applied across all samples.
#' @param min_abundance `Numeric` value. Minimum abundance threshold for a feature to be retained.
#'   Must be non-negative. Features with abundance below this are considered absent.
#' @param min_prevalence `Numeric` value. Minimum prevalence required for retention.
#'   If value is < 1, interpreted as proportion of samples; otherwise, as absolute number of samples.
#' @param abundance_criterion `Character` string. Criterion to use for filtering:
#'   \describe{
#'     \item{\code{prevalence}:}{Retain features present in at least \code{min_prevalence} samples
#'       (within group if \code{group_var} is used) and with abundance >= \code{min_abundance}
#'       in those samples.}
#'     \item{\code{mean}:}{Also requires that the mean abundance across samples (or group)
#'       is >= \code{min_abundance}.}
#'   }
#'   Default: \code{"prevalence"}.
#' @param group_requirement `Character` string. When `group_var` is specified, determines
#'   whether the criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#' @param keep_filtered `Logical`. If `TRUE` and any features are removed, stores the
#'   otu_table and tax_table of the filtered-out features in the `@info$filtered_taxa`
#'   slot of `microEDA` objects as matrices. Has no effect for `phyloseq` objects.
#'
#' @return A new `microEDA` or `phyloseq` object with features not meeting the
#'   requirements removed. A message reports how many features were dropped. If
#'   `keep_filtered = TRUE` and features were removed, the abundances and taxonomy
#'   of the removed features are stored in the `filtered_taxa()` slot of `microEDA`
#'   objects.
#'
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' me <- microEDA(GlobalPatterns)
#'
#' # Filter features present in at least 10% of samples with an abundance >= 100
#' filtered_me <- filter_features(me, min_abundance = 100, min_prevalence = 0.1)
#'
#' # Stratified filtering: feature must meet criteria in all groups, that is
#' # in at least 2 samples per group with an abundance of >= 10
#' filtered_me <- filter_features(
#'   me,
#'   group_var = "SampleType",
#'   min_abundance = 10,
#'   min_prevalence = 2,
#'   group_requirement = "all"
#' )
#'
#' @export
#' @importFrom phyloseq otu_table sample_data prune_taxa tax_table ntaxa
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr n inner_join group_by summarise mutate filter pull
#' @importFrom rlang .data
filter_features <- function(me,
                            min_abundance = 0,
                            min_prevalence = 0,
                            group_var = NULL,
                            abundance_criterion = c("prevalence", "mean"),
                            group_requirement = c("any", "all"),
                            keep_filtered = TRUE) {
  if (!inherits(me, "microEDA") && !inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }
  abundance_criterion <- match.arg(abundance_criterion, choices = c("prevalence", "mean"))
  group_requirement <- match.arg(group_requirement, choices = c("any", "all"))

  # Validate
  if (!is.numeric(min_abundance) || length(min_abundance) != 1 || min_abundance < 0) {
    stop("'min_abundance' must be a single non-negative numeric value")
  }

  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0) {
    stop("'min_prevalence' must be a single non-negative numeric value")
  }

  # Check for nonsensical filtering conditions
  if ((min_abundance == 0) != (min_prevalence == 0)) {
    warning("One of 'min_abundance' or 'min_prevalence' is zero while the other is not. This may lead to unintended filtering behavior as the zero threshold will pass all values.")
  }

  abund_data <- as(phyloseq::otu_table(me), "matrix")
  taxa_names <- rownames(abund_data)

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
      pass_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(pass_threshold, na.rm = TRUE)
      keep <- prevalence >= min_samples
    } else { # "mean"
      avg_abund <- rowMeans(abund_data, na.rm = TRUE)
      pass_threshold <- abund_data >= min_abundance
      prevalence <- rowSums(pass_threshold, na.rm = TRUE)
      keep <- avg_abund >= min_abundance & prevalence >= min_samples
    }
    taxa_to_keep <- taxa_names[keep]
  } else {
    # Grouped filtering
    metadata <- phyloseq::sample_data(me) |>
      tibble::rownames_to_column("Sample")

    if (!group_var %in% colnames(metadata)) stop("'group_var' not found in sample_data.")

    df_long <- tibble::rownames_to_column(as.data.frame(abund_data), "tax_ID") |>
      tidyr::pivot_longer(
        cols = -"tax_ID",
        names_to = "Sample", values_to = "Abundance"
      ) |>
      dplyr::inner_join(as.data.frame(metadata), by = "Sample")

    filter_fn <- if (group_requirement == "all") all else any

    if (abundance_criterion == "prevalence") {
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], .data$tax_ID) |>
        dplyr::summarise(
          pass = sum(.data$Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(.data$n_samples * min_prevalence))) |>
        dplyr::group_by(.data$tax_ID) |>
        dplyr::filter(filter_fn(.data$pass >= .data$min_n)) |>
        dplyr::pull(.data$tax_ID)
    } else { # "mean"
      keep <- df_long |>
        dplyr::group_by(.data[[group_var]], .data$tax_ID) |>
        dplyr::summarise(
          avg_abund = mean(.data$Abundance, na.rm = TRUE),
          pass = sum(.data$Abundance >= min_abundance),
          n_samples = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(min_n = ifelse(min_prevalence >= 1, as.integer(min_prevalence), ceiling(.data$n_samples * min_prevalence))) |>
        dplyr::group_by(.data$tax_ID) |>
        dplyr::filter(filter_fn(.data$avg_abund >= min_abundance & .data$pass >= .data$min_n)) |>
        dplyr::pull(.data$tax_ID)
    }
    taxa_to_keep <- unique(keep)
  }

  if (!length(taxa_to_keep) > 0) stop("All features filtered out!")

  filtered_data <- phyloseq::prune_taxa(taxa_to_keep, me)

  # Report on removed taxa
  n_removed <- phyloseq::ntaxa(me) - phyloseq::ntaxa(filtered_data)
  if (n_removed > 0) {
    message("Total filtered features: ", n_removed, "\n")
    # Extract filtered-out taxa
    dropped_taxa <- setdiff(taxa_names, taxa_to_keep)
    message("Filtered features: ", toString(dropped_taxa), "\n")
  }

  # Add "Other" row if taxa were removed
  if (keep_filtered && n_removed > 0) {
    if (inherits(me, "microEDA")) {
      # After filtering, store only minimal data to minimize object size
      filtered_otu <- otu_table(me)[dropped_taxa, , drop = FALSE]
      filtered_tax <- tax_table(me)[dropped_taxa, , drop = FALSE]

      # Retrieve existing filtered taxa, if any
      existing_data <- filtered_taxa(filtered_data)

      # Merge with existing data if present
      if (!is.null(existing_data)) {
        if (!is.null(existing_data$filtered_otu) && nrow(existing_data$filtered_otu) > 0) {
          filtered_otu <- rbind(existing_data$filtered_otu, filtered_otu)
        }
        if (!is.null(existing_data$filtered_tax) && nrow(existing_data$filtered_tax) > 0) {
          filtered_tax <- rbind(existing_data$filtered_tax, filtered_tax)
        }
      }

      # Save as matrices in microEDA@info, avoiding full phyloseq overhead
      filtered_taxa(filtered_data) <- list(
        filtered_otu = as(filtered_otu, "matrix"),
        filtered_tax = as(filtered_tax, "matrix")
      )
    }
  }

  if (inherits(me, "microEDA")) {
    # Reconstruct microEDA
    filtered_obj <- new("microEDA", filtered_data)
    if (n_removed > 0) {
      filter_history(filtered_obj) <- list(
        min_abundance = min_abundance,
        min_prevalence = min_prevalence,
        group_var = group_var,
        abundance_criterion = abundance_criterion,
        group_requirement = group_requirement,
        keep_filtered = keep_filtered,
        n_removed = n_removed
      )
    }
  } else {
    filtered_obj <- filtered_data # As it will already be a phyloseq object in this case
  }
  return(filtered_obj)
}


#' @param me A `microEDA` object containing filtered taxa in filtered_taxa(me).
#'
#' @returns A `list` with two elements: other_row, a 1-row matrix/otu_table of
#'  collapsed/summed abundances for filtered taxa merged into "Other", and other_tax,
#'  a 1-row matrix/tax_table assigning the "Other" label to all taxonomic ranks.
#'  Returns `NULL` if `me` is not a `microEDA` object or has no filtered taxa.
#'
#' @keywords internal
#' @noRd
.collapse_other <- function(me) {
  if (!inherits(me, "microEDA") || is.null(filtered_taxa(me))) {
    return(NULL)
  }

  # To reduce object size, avoid full phyloseq overhead and use matrices instead
  otu <- filtered_taxa(me)$filtered_otu
  tax <- filtered_taxa(me)$filtered_tax

  # Collapse abundances
  other_row <- matrix(colSums(otu), nrow = 1, dimnames = list("Other", colnames(otu)))
  other_tax <- matrix("Other",
    nrow = 1, ncol = ncol(tax),
    dimnames = list("Other", colnames(tax))
  )
  list(other_row = other_row, other_tax = other_tax)
}


#' @keywords internal
#' @noRd
.filter_threshold <- function(df, threshold, tax_col = "tax_ID") {
  # Identify numeric columns (sample/abundance columns)
  sample_cols <- colnames(df)[vapply(df, is.numeric, logical(1L)) & colnames(df) != tax_col]

  # Check if "Other" already exists
  other_row <- which(df[[tax_col]] == "Other")

  # Remove existing "Other" row if present, save it for reuse
  if (length(other_row) > 0) {
    other_data <- df[other_row, , drop = FALSE]
    df <- df[-other_row, ]
  } else {
    other_data <- data.frame(lapply(df, function(x) {
      if (is.numeric(x)) 0 else ""
    }), stringsAsFactors = FALSE)
    other_data[[tax_col]] <- "Other"
    rownames(other_data) <- "Other"
  }

  # Find rows below threshold (excluding any previous "Other")
  below_threshold <- apply(df[sample_cols], 1, max, na.rm = TRUE) < threshold

  # Sum abundances into the "Other" row, preserving existing values if present
  other_data[sample_cols] <- other_data[sample_cols] + colSums(df[below_threshold, sample_cols], na.rm = TRUE)

  # Remove rows below threshold
  df <- df[!below_threshold, ]

  # Append "Other" row only if it has non-zero values
  if (any(other_data[sample_cols] > 0)) {
    # df <- rbind(df, other_data)
    df <- dplyr::bind_rows(df, other_data)
  }

  return(df)
}


#' Apply Incrementing Relative Abundance Filter
#'
#' Iteratively increases a relative abundance threshold until no more than `ntaxa`
#' taxa remain above the threshold. Taxa below the threshold are merged into an
#' "Other" category. Works on both count and relative abundance data.
#'
#' @param df Data frame with taxonomic IDs and numeric abundance (samples) columns.
#' @param ntaxa Maximum number of taxa to retain before merging into "Other".
#' @param tax_column Name of the column containing taxonomic IDs.
#' @param initial_threshold Starting threshold for relative abundance (%).
#' @param increment Amount to increase threshold per iteration.
#' @param verbose Print message with final threshold used.
#' @return List containing filtered data and final threshold.
#' @details
#' The function iteratively increases the relative abundance threshold until no
#' more than `ntaxa` unique taxa remain. Taxa falling below the threshold in
#' each iteration are aggregated into an "Other" category.
#' Note that the final abundance of "Other" in any sample may exceed the reported
#' threshold because it is the sum of multiple individual taxa, each of which had
#' a relative abundance below the threshold.
#' @keywords internal
#' @noRd
.apply_incrementing_filter <- function(df,
                                       ntaxa = 30,
                                       tax_column = "tax_ID",
                                       initial_threshold = 0,
                                       increment = 0.5,
                                       verbose = TRUE) {
  initial_data <- df

  # Identify sample (abundance) columns
  sample_cols <- names(df)[vapply(df, is.numeric, logical(1L)) & names(df) != tax_column]

  # Convert to relative abundance (%)
  rel_data <- df
  rel_data[sample_cols] <- proportions(as.matrix(df[sample_cols]), margin = 2) * 100

  # Count non-"Other" taxa
  n_unique_taxa <- length(setdiff(rel_data[[tax_column]], "Other"))

  # Iteratively increase threshold until <= ntaxa taxa remain
  threshold <- initial_threshold
  while (n_unique_taxa > ntaxa) {
    rel_data <- .filter_threshold(rel_data, threshold = threshold, tax_col = tax_column)
    # Update n_unique_taxa but exclude "Other"
    n_unique_taxa <- length(setdiff(rel_data[[tax_column]], "Other"))
    threshold <- threshold + increment
  }

  # Get list of taxa to keep in original data (excluding those merged into
  # "Other" if present)
  kept_taxa <- setdiff(rel_data[[tax_column]], "Other")

  # Subset original data to keep only these taxa
  kept_rows <- initial_data[initial_data[[tax_column]] %in% kept_taxa, ]

  # Only create "Other" if some taxa were excluded
  if (length(kept_taxa) < length(unique(initial_data[[tax_column]]))) {
    # Sum abundances of all excluded taxa per sample to create "Other"
    other_mask <- !initial_data[[tax_column]] %in% kept_taxa
    other_abund <- colSums(initial_data[other_mask, sample_cols, drop = FALSE], na.rm = TRUE)

    other_row <- data.frame(as.list(other_abund), check.names = FALSE, stringsAsFactors = FALSE)
    other_row[[tax_column]] <- "Other"
    rownames(other_row) <- "Other"

    # Combine kept taxa and "Other" row
    filtered_data <- dplyr::bind_rows(kept_rows, other_row)
  } else {
    filtered_data <- kept_rows
  }

  # Rename "Other" to reflect threshold
  other_index <- which(filtered_data[[tax_column]] == "Other")
  if (length(other_index) > 0) {
    filtered_data[[tax_column]][other_index] <- paste0("Other (<", threshold, "%)")
  }

  if (verbose && threshold > 0) {
    message("Taxa below ", threshold, "% relative abundance were merged into 'Other'.")
  }

  # Return subset data and final 'Other' threshold
  return(list(
    data = filtered_data,
    threshold = if (threshold > 0) threshold else NULL
  ))
}
