#' Check if all columns in a profile are integer
#'
#' @param tax_profile A data frame or matrix to check.
#' @return `TRUE` if all columns contain only integer values (ignoring `NA`), `FALSE` otherwise.
#' @keywords internal
#' @noRd
.check_if_int <- function(tax_profile) {
  if (!all(vapply(tax_profile, is.numeric, logical(1L)))) stop("Non-numeric columns detected.")
  all(vapply(tax_profile, function(col) {
    all(col == as.integer(col), na.rm = TRUE)
  }, logical(1L)))
}


#' Check if a taxonomic profile contains proportional data
#'
#' @param tax_profile A matrix or data frame containing taxonomic abundance data.
#' Non-numeric columns are ignored.
#' @param tolerance Numeric threshold for acceptable variation in column sums.
#' @param silent `Logical`; if FALSE (default), warnings are issued for failed checks.
#' If TRUE, the function runs silently.
#'
#' @details
#' Determines whether the numeric columns in a taxonomic profile represent
#' relative abundances (proportions) by verifying that: (1) column sums are
#' approximately constant within a tolerance, (2) no negative values are present,
#' and (3) all values are within valid proportion ranges (between 0 and 1, or 0 and 100).
#'
#' @returns `TRUE` if data meet all criteria for relative abundances; `FALSE` otherwise.
#' @noRd
.is_proportion <- function(tax_profile, tolerance = 0.3, silent = FALSE) {
  is_rel <- TRUE

  is_numeric <- if (is.matrix(tax_profile)) {
    rep(is.numeric(tax_profile), ncol(tax_profile))
  } else {
    vapply(tax_profile, is.numeric, logical(1))
  }
  df_num <- tax_profile[, is_numeric, drop = FALSE]

  col_sums <- colSums(df_num, na.rm = TRUE)
  is_constant <- diff(range(col_sums)) < tolerance

  if (!is_constant) {
    if (!silent) warning("Column sums vary - data may not be relative abundances or may be transformed.")
    is_rel <- FALSE
  }

  # Check for negative values
  if (any(df_num < 0, na.rm = TRUE)) {
    if (!silent) warning("Negative values present - data likely log-transformed.")
    is_rel <- FALSE
  }

  # Check if values are proportions (0 ≤ x ≤ 1(00))
  # Upper limit needs to be a bit higher to account for minor floating-point
  # precision errors
  if (!all(df_num >= 0 & df_num <= 1.02 | df_num <= 102, na.rm = TRUE)) {
    if (!silent) warning("Values outside [0,1] or [0,100] - data may not be relative abundances.")
    is_rel <- FALSE
  }
  return(is_rel)
}


#' Check if a taxonomic profile contains raw count data
#'
#' @param tax_profile A matrix or data frame containing taxonomic abundance data.
#' Non-numeric columns are ignored.
#' @param tolerance Numeric threshold for acceptable variation in column sums.
#' @param silent `Logical`; if FALSE (default), warnings are issued for failed checks.
#' If TRUE, the function runs silently.
#'
#' @details
#' Determines whether the numeric columns in a taxonomic profile represent
#' raw count data by verifying that: (1) values are non-negative integers,
#' (2) column sums vary significantly (not constant like proportions),
#' (3) values are outside typical proportion ranges (between 0 and 1, or 0 and 100),
#' and (4) no negative or decimal values indicative of transformation (e.g., log, CLR).
#'
#' @returns `TRUE` if data meet all criteria for raw count data; `FALSE` otherwise.
#' @noRd
.is_counts <- function(tax_profile, tolerance = 0.3, silent = FALSE) {
  is_count <- TRUE

  # Identify numeric columns
  is_numeric <- if (is.matrix(tax_profile)) {
    rep(is.numeric(tax_profile), ncol(tax_profile))
  } else {
    vapply(tax_profile, is.numeric, logical(1))
  }
  df_num <- tax_profile[, is_numeric, drop = FALSE]

  # Check for non-integer values (indicative of transformation)
  if (any(vapply(df_num, function(x) any(x != as.integer(x), na.rm = TRUE), logical(1L)), na.rm = TRUE)) {
    if (!silent) warning("Decimal values present - data may be proportions or transformed.")
    is_count <- FALSE
  }

  # Check for negative values
  if (any(df_num < 0, na.rm = TRUE)) {
    if (!silent) warning("Negative values present - data is likely transformed.")
    is_count <- FALSE
  }

  # Check if values are within proportion ranges (0–1 or 0–100)
  if (all(df_num >= 0 & df_num <= 1.02, na.rm = TRUE) || all(df_num >= 0 & df_num <= 102, na.rm = TRUE)) {
    if (!silent) warning("Values within [0,1] or [0,100] range - data may be proportions.")
    is_count <- FALSE
  }

  # Check if column sums are approximately constant (expected in proportions, not counts)
  col_sums <- colSums(df_num, na.rm = TRUE)
  is_constant <- diff(range(col_sums)) < tolerance
  if (is_constant) {
    if (!silent) warning("Column sums are nearly constant - data may be proportions or transformed.")
    is_count <- FALSE
  }
  return(is_count)
}


#' Check if object is a valid taxonomic profile
#'
#' @param tax_profile A data frame with taxonomic abundance values for sample data.
#' @param tax_col Character vector of column names expected to contain taxonomic information (e.g. "clade_name").
#' @return `TRUE` if at least one taxonomic column exists and all others are numeric, `FALSE` otherwise.
#' @keywords internal
#' @noRd
.check_if_profile <- function(tax_profile, tax_col) {
  # Get columns with taxonomic information
  # taxcol_id <- unlist(lapply(tax_col, grep, colnames(tax_profile)))
  taxcol_id <- which(names(tax_profile) %in% tax_col) # actually a bit faster

  # Get sample columns
  sample_cols <- tax_profile[, -taxcol_id, drop = FALSE]

  # Check that there is a taxonomic column and check that rest of the columns
  # represent abundances in samples.
  isTRUE(length(taxcol_id) > 0 &&
    all(vapply(sample_cols, is.numeric, logical(1L))))
}


#' Check if Taxonomic Rank is Valid
#'
#' Internal function to validate whether a given taxonomic level is one of the
#' accepted ranks. The function checks against a predefined set of valid ranks,
#' including both abbreviated (e.g., "k", "p") and full names (e.g., "kingdom", "Phylum").
#'
#' @param taxa_lvl A `character` string representing the taxonomic rank to be validated.
#'
#' @return A `logical` value: `TRUE` if the taxonomic rank is valid, `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' .is_valid_rank("Genus") # TRUE
#' .is_valid_rank("g") # TRUE
#' .is_valid_rank("tribe") # FALSE
.is_valid_rank <- function(taxa_lvl) {
  stopifnot(is.character(taxa_lvl))

  taxa_lvl <- tolower(taxa_lvl)

  valid_taxa_lvls <- c(
    "k", "p", "c", "o", "f", "g", "s", "t", "kingdom", "phylum",
    "class", "order", "family", "genus", "species", "strain"
  )
  taxa_lvl %in% valid_taxa_lvls
}


#' Valid Taxonomic Ranks Error Message
#'
#' Internal character vector containing the error message that should be returned
#' when an invalid taxonomic rank is provided. Intended to be used in conjunction
#' with `.is_valid_rank()` for input validation in other functions.
#'
#' @keywords internal
#' @noRd
.valid_ranks_msg <- "Invalid tax_rank. Please choose one of 'kingdom', 'phylum','class', 'order', 'family', 'genus', 'species' or 'strain'."


#' Check for reserved column names and rename if necessary
#'
#' @param x The data frame to check.
#' @param reserved_var_names `Character` vector of reserved names (default: c("Sample", "Abundance", "tax_ID")).
#' @return The input data frame with renamed columns if conflicts exist, and a warning.
#' @keywords internal
#' @noRd
.check_var_names <- function(x, reserved_var_names = c("Sample", "Abundance", "tax_ID")) {
  vars <- colnames(x)

  var_conflict <- intersect(reserved_var_names, vars)

  if (length(var_conflict) > 0) {
    idx <- which(vars %in% var_conflict)
    new_names <- paste0("orig_", vars[idx])
    msg <- paste(
      "The variables: \n",
      paste(vars[idx], collapse = ", "),
      "\n  have been renamed to:\n",
      paste0(new_names, collapse = ", "),
      "\n  to avoid conflicts with special plot attribute names."
    )
    .show_warning(msg)
    colnames(x)[idx] <- new_names
  }
  return(x)
}


#' Display warning in call stack
#'
#' @description This function displays the given `warnmessage` as a warning
#'  message. Optionally users can specify how many frames back in the call
#'  sequence this should be shown (default: top-level function that initiated
#'  the call sequence).
#'
#' @param warnmessage The warning message to display. Needs to be one string.
#' @param n `Integer` value indicating how many frames back in the call sequence
#'  the warning should be shown. Default is NULL, showing the warning for the
#'  top-level function that initiated the call sequence.
#'
#' @returns The warning message as `character` string, invisibly.
#' @keywords internal
#' @noRd
#'
#' @examples
#' a <- function() b()
#' b <- function() {
#'   c()
#'   d()
#' }
#' c <- function() .show_warning("This will display a warning for the top-level function")
#' d <- function() .show_warning("This will display a warning for the 'b' function", n = 2)
#' a()
.show_warning <- function(warnmessage, n = NULL) {
  if (is.numeric(n)) {
    calling_frame <- sys.call(sys.nframe() - n)
  } else {
    # This returns the name of the top-level function that initiated the call sequence
    calling_frame <- sys.calls()[[1]]
  }
  call_str <- paste(deparse(calling_frame), collapse = " ")
  # Reduce multiple spaces to one for long messages due to deparse() preserving
  # spacing from the parse tree
  call_str <- gsub("\\s+", " ", call_str)
  # msg <- paste("In", paste(deparse(calling_frame), collapse = " "), ":\n ", warnmessage)
  msg <- paste("In", call_str, ":\n ", warnmessage)
  warning(msg, call. = FALSE)
}


#' Display error message in call stack
#'
#' @description his function displays the given `errormessage` as a warning
#'  message. Optionally users can specify how many frames back in the call
#'  sequence this should be shown (default: top-level function that initiated
#'  the call sequence).
#'
#' @param errormessage The error message to display. Must be a single string.
#' @param n `Integer` indicating how many frames back in the call sequence
#'  the error should be shown. Default is NULL, showing the error for the
#'  top-level function that initiated the call sequence.
#'
#' @returns (Invisibly) the error message as a `character` string, though
#'  execution typically halts due to `stop()`.
#' @keywords internal
#' @noRd
#'
#' @examples
#' a <- function() b()
#' b <- function() {
#'   c()
#'   d()
#' }
#' c <- function() .show_error("Error at top level")
#' d <- function() .show_error("Error attributed to 'b'", n = 2)
#' # a()  # Uncomment to test
.show_error <- function(errormessage, n = NULL) {
  if (is.numeric(n)) {
    calling_frame <- sys.call(sys.nframe() - n)
  } else {
    # This returns the name of the top-level function that initiated the call sequence
    calling_frame <- sys.calls()[[1]]
  }
  call_str <- paste(deparse(calling_frame), collapse = " ")
  # Reduce multiple spaces to one for long messages due to deparse() preserving
  # spacing from the parse tree
  call_str <- gsub("\\s+", " ", call_str)
  # msg <- paste("In", paste(deparse(calling_frame), collapse = " "), ":\n ", warnmessage)
  msg <- paste("In", call_str, ":\n ", errormessage)
  stop(msg, call. = FALSE)
}


#' Warn about and filter disallowed arguments in ellipsis
#'
#' @description This function is used to restrict which arguments can be set
#'  through the ellipsis (`...`). It checks for allowed parameters and warns
#'  on disallowed ones while removing them from the `arglist` list.
#'
#' @param allowed A `character` vector specifying the names of allowed parameters.
#' @param arglist A named `list` containing all provided arguments in the ellipsis (`...`).
#'
#' @returns The modified `arglist` list with only allowed arguments.
#' @keywords internal
#' @noRd
#'
#' @examples
#' example_fun <- function(...) {
#'   allowed_args <- c("fill", "alpha")
#'   arglist <- list(...)
#'   .warn_invalid_args(allowed_args, arglist)
#' }
#' example_fun(fill = "red", color = "blue", alpha = 0.7)
.warn_invalid_args <- function(allowed, arglist) {
  provided <- names(arglist)

  # Check if any arguments are not in the allowed vector
  disallowed <- setdiff(provided, allowed)
  if (length(disallowed)) {
    .show_warning(warnmessage = paste(
      "Ignoring invalid arguments:",
      paste(disallowed, collapse = ", ")
    ))
  }
  # Keep only allowed args
  arglist <- arglist[intersect(provided, allowed)]
}


#' Check if sample (meta)data is present in object
#'
#' @param x A `microEDA` or `phyloseq` object
#'
#' @returns `TRUE` if sample data is present; `FALSE` otherwise.
#' @keywords internal
#' @noRd
.check_sample_data <- function(x) {
  stopifnot(inherits(x, "phyloseq"))
  metadata <- phyloseq::sample_data(x, errorIfNULL = FALSE)
  # If no metadata the above will simply return NULL
  length(metadata) > 0
}


#' Check Taxonomic Consistency Across Higher-Level Classifications
#'
#' Verifies that each taxon at a specified rank has a consistent higher-level
#' taxonomic path (e.g., all features labeled as the same genus belong to the
#' same family, order, etc.).
#' Useful for identifying inconsistencies (from, e.g., different database versions)
#' in taxonomic annotations within microbiome datasets.
#'
#' @param me A `microEDA` or `phyloseq` object, or a `data.frame/matrix` where rows
#'   represent taxa and columns represent taxonomic ranks (e.g., Kingdom, Phylum,
#'   Class, etc.).
#' @param tax_rank `Character` string specifying the taxonomic rank at which to
#'   check consistency (e.g., "Genus", "Family"). If `NULL`, defaults to the
#'   lowest (last) taxonomic rank present in the data.
#' @param detailed_report `Logical`. If `TRUE`, provides a detailed report showing
#'   at which higher taxonomic levels inconsistencies occur. If `FALSE`, only
#'   returns the names of inconsistent taxa.
#'
#' @return Invisibly returns one of the following:
#'   \itemize{
#'     \item `NULL` if no inconsistencies are found.
#'     \item A `character` vector of inconsistent taxon names if `detailed_report = FALSE`.
#'     \item A `tibble` with detailed per-rank inconsistency flags if
#'           `detailed_report = TRUE` and inconsistencies exist.
#'   }
#'   To use the return value, assign the function call to a variable (e.g.,
#'   \code{result <- check_taxonomic_consistency(me, detailed_report = TRUE)}).
#'
#' @details
#' The function constructs a "clade lineage" by concatenating all taxonomic ranks up to
#' the specified `tax_rank`. It then groups by the `tax_rank` value and checks whether
#' multiple clade paths exist for the same taxon at that rank. If so, a warning is issued.
#'
#' When `detailed_report = TRUE`, the function performs additional analysis to
#' identify exactly which higher ranks differ for each inconsistent taxon.
#'
#' @examples
#' # Example using a simple data frame
#' tax_data <- data.frame(
#'   Kingdom = rep("Bacteria", 4),
#'   Phylum = c("Firmicutes", "Firmicutes", "Proteobacteria", "Firmicutes"),
#'   Class = c("Bacilli", "Bacilli", "Gammaproteobacteria", "Bacilli"),
#'   Genus = c("Lactobacillus", "Lactobacillus", "Escherichia", "Lactobacillus")
#' )
#'
#' # Check consistency at Genus level
#' check_taxonomic_consistency(tax_data, tax_rank = "Genus")
#' @importFrom rlang .data
#' @export
check_taxonomic_consistency <- function(me, tax_rank = NULL, detailed_report = TRUE) {
  # Input validation
  if (inherits(me, "microEDA") || inherits(me, "phyloseq")) {
    tax_df <- as.data.frame(phyloseq::tax_table(me))
  } else if (is.data.frame(me) || is.matrix(me)) {
    tax_df <- as.data.frame(me)
  } else {
    stop("'me' must be a microEDA, phyloseq, data.frame, or matrix object.")
  }

  # Default to lowest (last) taxonomic rank if not specified
  if (is.null(tax_rank)) {
    rank_names <- colnames(tax_df)
    tax_rank <- rank_names[length(rank_names)]
  } else {
    tax_rank <- .get_full_tax_rank(tax_rank) # In case it was abbreviated
  }

  # Ensure rank exists
  if (!(tax_rank %in% colnames(tax_df))) {
    stop("Taxonomic rank '", tax_rank, "' not found in tax_table.")
  }

  ranks_up_to <- names(tax_df)[1:which(names(tax_df) == tax_rank)]

  # Create full clade path
  tax_df <- tidyr::unite(
    tax_df,
    "clade_name",
    dplyr::all_of(ranks_up_to),
    sep = ";",
    remove = FALSE
  )

  # Group by target rank and check for multiple clade paths
  inconsistencies <- tax_df |>
    dplyr::group_by(.data[[tax_rank]]) |>
    dplyr::mutate(all_same = dplyr::n_distinct(.data$clade_name) == 1) |>
    dplyr::filter(!.data$all_same) |>
    dplyr::distinct(.data$clade_name, .keep_all = TRUE) |>
    dplyr::ungroup()

  if (nrow(inconsistencies) == 0) {
    message("No taxonomic inconsistencies detected at rank '", tax_rank, "'.")
    return(invisible(NULL))
  }

  inconsistent_taxa <- dplyr::pull(inconsistencies, .data[[tax_rank]]) |> unique()
  n_inconsistent <- length(inconsistent_taxa)
  taxon_label <- if (n_inconsistent == 1) "taxon" else "taxa"

  # Temporarily enable immediate warning print so it appears before detailed_report
  old_warn <- options(warn = 1)
  on.exit(options(old_warn), add = TRUE)

  warning(
    "Conflicting taxonomy for ", n_inconsistent, " ", taxon_label,
    " at rank '", tax_rank, "'. Inconsistent higher-level classification",
    if (n_inconsistent > 1) "s", " detected for:\n",
    toString(inconsistent_taxa),
    call. = FALSE
  )

  # Enhanced diagnostic: see at which ranks taxonomic inconsistencies were found
  # for a taxon
  if (detailed_report && nrow(inconsistencies) > 0) {
    detailed_differences <- inconsistencies |>
      dplyr::group_by(.data[[tax_rank]]) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::any_of(ranks_up_to),
          ~ dplyr::n_distinct(.x) > 1,
          .names = "differs_in_{col}"
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(dplyr::if_any(dplyr::starts_with("differs_in_"), ~ .x == TRUE))

    # Print detailed differences
    if (nrow(detailed_differences) > 0) {
      message("Detailed differences (TRUE = inconsistent rank):")
      print(detailed_differences)
      # Return detailed info for further inspection if too many rows
      return(invisible(detailed_differences))
    }
  }
  # Otherwise return just the names of inconsistent taxa
  return(invisible(inconsistent_taxa))
}
