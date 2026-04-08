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
