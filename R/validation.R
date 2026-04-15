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
    if (!silent) warning("Values outside [0,1] or [0,100] - may not be relative abundances.")
    is_rel <- FALSE
  }
  return(is_rel)
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
