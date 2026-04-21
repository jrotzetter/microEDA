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


#' Filter Features in a microEDA Object Based on Abundance and Prevalence
#'
#' Filters taxa (features) from a `microEDA` object based on minimum abundance
#' and prevalence thresholds. Filtering can be applied globally or within groups
#' defined by a metadata variable. Optionally adds an "Other" category summarizing
#' filtered features.
#'
#' @param me A `microEDA` object containing microbiome abundance data and associated metadata.
#' @param group_var (optional) Name of the sample metadata variable to define groups
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
#' @param keep_other `Logical`. If `TRUE` and any features are removed, adds a new row
#'   named `"Other"` to the OTU and tax tables, summarizing the total abundance of
#'   all filtered features. Existing `"Other"` rows are merged.
#'
#' @return A new `microEDA` object with filtered features. A message reports how many
#'   features were removed. If `keep_other = TRUE` and features were removed,
#'   an `"Other"` row is added or updated.
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
filter_features <- function(me,
                            min_abundance = 0,
                            min_prevalence = 0,
                            group_var = NULL,
                            abundance_criterion = c("prevalence", "mean"),
                            group_requirement = c("any", "all"),
                            keep_other = TRUE) {
  if (!inherits(me, "microEDA")) {
    stop("'me' must be a microEDA object")
  }
  abundance_criterion <- match.arg(abundance_criterion)
  group_requirement <- match.arg(group_requirement)

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

    df_long <- tibble::rownames_to_column(as.data.frame(abund_data), "tax_ID") |>
      tidyr::pivot_longer(
        cols = -"tax_ID",
        names_to = "Sample", values_to = "Abundance"
      ) |>
      dplyr::inner_join(as.data.frame(metadata), by = "Sample")

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
    taxa_to_keep <- unique(keep)
  }

  if (!length(taxa_to_keep) > 0) stop("All features filtered out!")

  filtered_data <- phyloseq::prune_taxa(taxa_to_keep, me)

  # Report on removed taxa
  n_removed <- phyloseq::ntaxa(me) - phyloseq::ntaxa(filtered_data)
  if (n_removed > 0) {
    message("Total filtered features: ", n_removed, "\n")
    dropped <- setdiff(taxa_names, taxa_to_keep)
    message("Filtered features: ", toString(dropped), "\n")
  }

  # Add "Other" row if taxa were removed
  if (keep_other && n_removed > 0) {
    orig_otu <- phyloseq::otu_table(me)
    pruned_otu <- phyloseq::otu_table(filtered_data)
    pruned_tax <- phyloseq::tax_table(filtered_data)

    # Check for existing 'Other' rows
    other_mask <- rownames(pruned_otu) == "Other"

    # Sum all 'Other' rows if any exist
    if (any(other_mask)) {
      other_sum <- colSums(pruned_otu[other_mask, , drop = FALSE], na.rm = TRUE)
      # Remove all existing 'Other' rows
      pruned_otu <- pruned_otu[!other_mask, ]
      pruned_tax <- pruned_tax[!other_mask, ]
    } else {
      other_sum <- 0
    }

    # Compute "Other" as difference
    other_vals <- colSums(orig_otu) - colSums(pruned_otu) + other_sum
    other_row <- matrix(other_vals,
      nrow = 1,
      dimnames = list("Other", colnames(other_vals))
    )

    # Extend OTU table
    new_otu <- rbind(pruned_otu, other_row)
    filtered_data@otu_table <- phyloseq::otu_table(new_otu, taxa_are_rows = TRUE)

    # Extend tax_table with "Other" for all taxonomic ranks
    n_tax_ranks <- ncol(pruned_tax)
    other_tax <- matrix("Other",
      nrow = 1, ncol = n_tax_ranks,
      dimnames = list("Other", colnames(pruned_tax))
    )

    # Bind to existing tax table
    new_tax <- rbind(pruned_tax, other_tax)
    filtered_data@tax_table <- phyloseq::tax_table(new_tax)
  }
  # Reconstruct microEDA
  microEDA_filtered <- new("microEDA", filtered_data, info = me@info)
}
