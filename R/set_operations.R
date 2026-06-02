#' Get a List of Present Taxa by Group
#'
#' Filters a `microEDA` or `phyloseq` object based on abundance and prevalence
#' criteria, then returns a named list of unique taxa present in each group
#' defined by a sample metadata variable. This function is useful for identifying
#' core microbiomes or group-specific taxa sets.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU table, taxonomy, and sample data.
#' @param tax_rank `Character` string specifying the taxonomic rank to analyze
#'                 (e.g., "Phylum", "Family"). Must be a valid rank in the
#'                 taxonomy table.
#' @param group_var `Character` string indicating a sample metadata
#'                  variable to group samples.
#' @param group_labels Named character `vector` mapping old group names to
#'                 new labels (e.g., `c("Old" = "New")`).
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
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param group_requirement `Character` string. When `group_var` is specified, determines
#'   whether the criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#' @param keep_filtered `Logical`. If `TRUE` and any features are removed, stores the
#'   otu_table and tax_table of the filtered-out features in the `@info$filtered_taxa`
#'   slot of `microEDA` objects as matrices. Has no effect for `phyloseq` objects.
#' @param rm_missing `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae").
#' @param add_prefix `Logical`. If `TRUE`, adds QIIME-style prefixes (e.g., `k__`, `p__`)
#'   to taxonomic labels.
#'
#' @return A named `list` where each element corresponds to a level of `group_var`.
#'   Each element contains a character vector of unique taxon names that meet the
#'   specified abundance and prevalence criteria.
#'
#' @details
#'   The function first validates the input `me` object and ensures the
#'   OTU table contains counts or proportions. It then processes the taxonomic
#'   profile, filtering taxa by `min_abundance`, `min_prevalence` using [filter_features],
#'   before agglomerating to the specified `tax_rank` with [agglomerate_taxa].
#'   Finally, it splits the resulting taxa by the grouping variable and returns
#'   a list of unique taxa per group.
#'
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' present <- get_presence_list(GlobalPatterns, "Phylum", "SampleType",
#'   min_abundance = 100, min_prevalence = 0.65
#' )
#'
#' @importFrom rlang sym
#' @importFrom dplyr arrange left_join rename group_by filter
#' @importFrom tibble rownames_to_column
#' @export
get_presence_list <- function(me,
                              tax_rank,
                              group_var,
                              group_labels = NULL,
                              min_abundance = 0,
                              min_prevalence = 0,
                              abundance_criterion = c("prevalence", "mean"),
                              filter_by_group = FALSE,
                              group_requirement = c("any", "all"),
                              keep_filtered = TRUE,
                              rm_missing = FALSE,
                              add_prefix = FALSE) {
  if (!inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }

  if (!.is_counts(otu_table(me), silent = TRUE) && !.is_proportion(otu_table(me), silent = TRUE)) {
    stop("'otu_table' is neither counts nor relative abundance - data was likely transformed.")
  }

  if (!.is_valid_rank(tax_rank)) {
    stop(.valid_ranks_msg)
  }
  tax_rank <- .get_full_tax_rank(tax_rank) # In case it was abbreviated

  if (.check_sample_data(me)) {
    if (!(group_var %in% names(phyloseq::sample_data(me)))) {
      stop("'group_var' not found in sample metadata.")
    }

    if (!is.null(group_labels)) {
      phyloseq::sample_data(me) <- .rename_values(phyloseq::sample_data(me), group_var, group_labels)
    }

    metadata <- data.frame(phyloseq::sample_data(me),
      stringsAsFactors = FALSE
    )[, group_var, drop = FALSE] |>
      .check_var_names() |>
      tibble::rownames_to_column("Sample")

    if (is.factor(metadata[[group_var]])) {
      metadata[[group_var]] <- as.character(metadata[[group_var]])
    }
    group_var_sym <- rlang::sym(group_var)
  } else {
    stop("Can't use 'group_var' without sample metadata!")
  }

  tax_abund <- .prepare_tax_profile(me,
    tax_rank = tax_rank,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    ntaxa = phyloseq::ntaxa(me),
    group_var = group_var,
    abundance_criterion = abundance_criterion,
    filter_by_group = filter_by_group,
    group_requirement = group_requirement,
    keep_filtered = keep_filtered,
    rm_missing = rm_missing,
    transform = "None",
    add_prefix = add_prefix,
    process_taxon = FALSE,
    calculate_prevalence = FALSE
  )

  min_abundance <- ifelse(min_abundance == 0, .Machine$double.eps, min_abundance)

  tax_presence <- tax_abund |>
    dplyr::arrange(.data$Sample) |>
    dplyr::left_join(metadata, by = "Sample") |>
    dplyr::rename(Group = .data[[group_var]]) |>
    dplyr::group_by(.data$Group, .data$Sample) |>
    dplyr::filter(.data$Abundance >= min_abundance) # |>
  # dplyr::group_by(.data$Group) |>
  # dplyr::reframe("Taxon" = unique(.data$Taxon)) |>
  # dplyr::mutate(Presence = TRUE)

  # Create a list to store taxon for each Variable Of Interest
  # sets_list <- tax_presence |>
  #   group_by(.data$Group) |>
  #   summarize(taxon_list = list(.data$Taxon)) |>
  #   ungroup()

  # Split by Group and extract unique taxa
  sets_list <- split(tax_presence$Taxon, tax_presence$Group)
  sets_list <- lapply(sets_list, unique)

  # Convert the list into a named list
  # sets_list <- stats::setNames(sets_list$taxon_list, sets_list$Group)
  sets_list <- stats::setNames(sets_list, names(sets_list))

  return(sets_list)
}


#' @keywords internal
#' @noRd
.get_specific_elements <- function(sets_list, sets1, sets2 = NULL, unite_sets1 = TRUE) {
  if (!is.list(sets_list)) stop("'sets_list' must be a list of character vectors")
  if (any(!vapply(sets_list, is.character, logical(1L)))) stop("All elements of 'sets_list' must be character vectors")

  # In case indices were used instead of names
  if (is.numeric(sets1)) sets1 <- names(sets_list)[sets1]
  if (is.numeric(sets2)) sets2 <- names(sets_list)[sets2]

  # Helper to handle edge cases and safely reduce with union or intersect
  safe_reduce <- function(lst, fun) {
    if (length(lst) == 0) {
      return(character(0L))
    }
    if (length(lst) == 1) {
      return(unlist(lst, use.names = FALSE))
    }
    purrr::reduce(lst, fun)
  }

  if (is.null(sets2)) {
    # Get the remaining sets
    sets2 <- names(sets_list)[!names(sets_list) %in% sets1]
  }

  set1 <- safe_reduce(sets_list[sets1], if (unite_sets1) base::union else base::intersect)
  set2 <- safe_reduce(sets_list[sets2], base::union)

  base::setdiff(set1, set2)
}


#' @keywords internal
#' @noRd
.get_intersection <- function(sets_list, sets = NULL) {
  if (!is.list(sets_list)) stop("'sets_list' must be a list of character vectors")
  if (any(!vapply(sets_list, is.character, logical(1L)))) stop("All elements of 'sets_list' must be character vectors")

  # if (is.null(sets)) {
  #   overlap <- purrr::reduce(sets_list, function(x, y) base::intersect(x, y))
  # } else {
  #   sets_slice <- sets_list[sets]
  #   overlap <- purrr::reduce(sets_slice, function(x, y) base::intersect(x, y))
  # }
  # return(overlap)
  input_sets <- if (is.null(sets)) sets_list else sets_list[sets]

  # Handle edge cases
  if (length(input_sets) == 0) {
    return(character(0))
  }
  if (length(input_sets) == 1) {
    return(input_sets[[1]])
  }

  purrr::reduce(input_sets, base::intersect)
}


#' @keywords internal
#' @noRd
.get_union <- function(sets_list, sets = NULL) {
  if (!is.list(sets_list)) stop("'sets_list' must be a list of character vectors")
  if (any(!vapply(sets_list, is.character, logical(1L)))) stop("All elements of 'sets_list' must be character vectors")

  # if (is.null(sets)) {
  #   combined_sets <- purrr::reduce(sets_list, function(x, y) base::union(x, y))
  # } else {
  #   sets_slice <- sets_list[sets]
  #   combined_sets <- purrr::reduce(sets_slice, function(x, y) base::union(x, y))
  # }
  # return(combined_sets)
  input_sets <- if (is.null(sets)) sets_list else sets_list[sets]

  # Handle edge cases
  if (length(input_sets) == 0) {
    return(character(0))
  }
  if (length(input_sets) == 1) {
    return(input_sets[[1]])
  }

  purrr::reduce(input_sets, base::union)
}


#' Compute Taxa Overlaps and Specific Sets Across Groups
#'
#' Filters a `microEDA` or `phyloseq` object to identify present taxa per group,
#' then calculates all possible intersections and group-specific taxa sets.
#' This function is useful for identifying core vs. unique microbiome components
#' across conditions and can serve as an alternative to an UpSet plot.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU table, taxonomy, and sample data.
#' @param tax_rank `Character` string specifying the taxonomic rank to analyze
#'                 (e.g., "Phylum", "Family"). Must be a valid rank in the
#'                 taxonomy table.
#' @param group_var `Character` string indicating a sample metadata variable to group samples.
#' @param group_labels Named character `vector` mapping old group names to
#'                 new labels (e.g., `c("Old" = "New")`).
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
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param group_requirement `Character` string. When `group_var` is specified, determines
#'   whether the criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#' @param keep_filtered `Logical`. If `TRUE` and any features are removed, stores the
#'   otu_table and tax_table of the filtered-out features in the `@info$filtered_taxa`
#'   slot of `microEDA` objects. Has no effect for `phyloseq` objects.
#' @param rm_missing `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor.
#' @param add_prefix `Logical`. If `TRUE`, adds QIIME-style prefixes (e.g., `k__`, `p__`)
#'   to taxonomic labels.
#'
#' @return A named `list` of character vectors where each element corresponds to a specific overlap pattern.
#'   Empty intersection results are omitted from the final list.
#'
#' @details
#'   This function computes taxa overlap patterns across sample groups in a microbiome dataset.
#'   It first calls [get_presence_list] to obtain filtered taxa sets per group based on abundance
#'   and prevalence criteria. It then calculates all possible set intersections and specific taxa.
#'
#'   \itemize{
#'     \item **Single group**: Returns only the taxa specific to that group (i.e., present in the
#'           group and absent in all others - if others exist). If only one group is provided overall,
#'           it simply returns that group's filtered taxa.
#'     \item **Multiple groups**: Iterates through all combination sizes from 1 to N (where N is the
#'           number of groups). For each combination:
#'           \itemize{
#'             \item If the combination includes **all groups**, the result is the **core microbiome**
#'                   (taxa present in every group), labeled `"all-intersection"`.
#'             \item If the combination includes **some but not all groups**, the result includes taxa
#'                   present in **all of the selected groups** and **absent from all others**. These are
#'                   labeled as `"<group1>;...;<groupk>-intersection"`.
#'             \item If the combination size is **1**, the result includes taxa **unique to that group**
#'                   (i.e., not found in any other group), labeled as `"<group>-specific"`.
#'           }
#'   }
#'
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' # Calculate overlaps for Phyla across SampleTypes
#' overlaps <- get_taxa_overlaps(GlobalPatterns, "Phylum", "SampleType",
#'   min_abundance = 100, min_prevalence = 3
#' )
#'
#' # View taxa specific to "Soil"
#' overlaps[["Soil-specific"]]
#'
#' # View taxa shared by all groups
#' overlaps[["all-intersection"]]
#'
#' @seealso [get_presence_list] for the underlying presence filtering logic.
#' @importFrom utils combn
#' @export
get_taxa_overlaps <- function(me,
                              tax_rank,
                              group_var,
                              group_labels = NULL,
                              min_abundance = 0,
                              min_prevalence = 0,
                              abundance_criterion = c("prevalence", "mean"),
                              filter_by_group = FALSE,
                              group_requirement = c("any", "all"),
                              keep_filtered = TRUE,
                              rm_missing = FALSE,
                              add_prefix = FALSE) {
  if (!inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }

  sets <- get_presence_list(me,
    tax_rank,
    group_var,
    group_labels = group_labels,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    abundance_criterion = abundance_criterion,
    filter_by_group = filter_by_group,
    group_requirement = group_requirement,
    keep_filtered = keep_filtered,
    rm_missing = rm_missing,
    add_prefix = add_prefix
  )

  if (length(sets) == 0) {
    warning("No intersections found with the specified criteria; returning an empty list.")
    return(list())
  }
  if (length(sets) == 1) {
    return(sets)
  }

  overlaps <- list()
  n <- length(sets)
  set_names <- names(sets)

  # Generate all combinations from 1 to n
  for (i in 1:n) {
    combs <- utils::combn(set_names, i, simplify = FALSE)
    for (comb in combs) {
      # For i == 1: specific to that set (intersection of one, minus all others)
      # For 1 < i < n: intersection of `comb`, not in any other set
      # For i == n: full intersection
      if (i == n) {
        result <- .get_intersection(sets, comb)
        key <- "all-intersection"
      } else {
        others <- set_names[!set_names %in% comb]
        result <- .get_specific_elements(sets, comb, others, unite_sets1 = FALSE)
        key <- paste0(paste0(comb, collapse = ";"), "-intersection")
        if (i == 1) key <- paste0(comb, "-specific")
      }
      if (length(result) > 0) {
        overlaps[[key]] <- result
      }
    }
  }
  return(overlaps)
}
