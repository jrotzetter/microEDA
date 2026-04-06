#' A class to represent a MetaPhlAn profile.
#'
#' This class stores database version information and taxonomic profile data from MetaPhlAn output.
#'
#' @slot mpa_version `Character` string indicating the MetaPhlAn database version used.
#' @slot mpa_taxProfile A `data.frame` containing the taxonomic profile (relative abundances of clades).
#' @export
setClass("metaphlanProfile",
  slots = list(mpa_version = "character", mpa_taxProfile = "data.frame"),
  prototype = list(
    mpa_version = character(),
    mpa_taxProfile = data.frame()
  )
)


#' Load MetaPhlAn profile
#'
#' @description
#' Function to load a MetaPhlAn profile from a file path.
#'
#' @param file_path A `character` string containing the path to the file to read.
#' @importFrom utils read.table
#' @importFrom rlang :=
#'
#' @returns An object of class `metaphlanProfile`.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "merged_abundance_table.txt", package = "microEDA")
#'
#' mtp <- load_metaphlan(file_path)
load_metaphlan <- function(file_path) {
  profile_type <- NULL

  lines <- tryCatch(
    readLines(file_path, n = 10),
    error = function(e) {
      stop(
        "Cannot read file: ", file_path,
        "\nPlease check that the file is a valid MetaPhlAn file"
      )
    }
  )

  # Extract the line containing the used database
  db_index <- startsWith(lines, "#mpa")
  db_version <- lines[db_index]

  if (!length(db_version) > 0) {
    warning("No MetaPhlAn database version found")
  }

  # Check if it is a single or merged metaphlan profile
  if (length(which(startsWith(lines, "#"))) > 1) {
    profile_type <- "single"
  } else {
    profile_type <- "merged"
  }

  switch(profile_type,
    "single" = {
      file_name <- tools::file_path_sans_ext(basename(file_path))

      mtprofile <- tryCatch(
        {
          read.table(file_path,
            header = FALSE,
            sep = "\t", # Needed when 'additional_species' column isn't empty
            comment.char = "#",
            check.names = FALSE
          )
        },
        error = function(condition) {
          stop(
            "Cannot read file: ", file_path,
            "\nPlease check that the file is a valid merged MetaPhlAn file"
          )
        }
      )

      # Find first line that starts not with a '#'
      start_row <- which(!startsWith(lines, "#"))[1]
      # The line before contains the commented out column names
      single_colnames <- lines[start_row - 1]
      single_colnames <- single_colnames |>
        sub(pattern = "#", replacement = "") |>
        strsplit(split = "\t") |>
        unlist()

      colnames(mtprofile) <- single_colnames

      if (length(colnames(mtprofile)) != 4 ||
        all(!c("clade_name", "relative_abundance") %in% colnames(mtprofile))) {
        stop(
          "Loaded file has wrong number of columns to be a single ",
          'MetaPhlAn profile!\nOnly the columns "clade_name",',
          '"NCBI_tax_id", "relative_abundance", "additional_species"',
          "should be present."
        )
      }
      mtprofile <- mtprofile[, c("clade_name", "relative_abundance")]

      mtprofile <- dplyr::rename(mtprofile, {{ file_name }} := "relative_abundance")

      result <- new("metaphlanProfile",
        mpa_version = db_version,
        mpa_taxProfile = mtprofile
      )
    },
    "merged" = {
      mtprofile <- tryCatch(
        {
          read.table(file_path,
            header = TRUE,
            comment.char = "#",
            check.names = FALSE
          )
        },
        error = function(condition) {
          stop(
            "Cannot read file: ", file_path,
            "\nPlease check that the file is a valid merged MetaPhlAn file"
          )
        }
      )

      # ID in Metaphlan v2, clade_name in Metaphlan > v2
      if (!.check_if_profile(mtprofile, c("clade_name", "ID"))) {
        stop(
          "Cannot read file: ", file_path,
          "\nPlease check that the file is a valid merged MetaPhlAn file"
        )
      }

      if ("ID" %in% names(mtprofile) && !"clade_name" %in% names(mtprofile)) {
        mtprofile <- dplyr::rename(mtprofile, "clade_name" = "ID")
      }

      result <- new("metaphlanProfile",
        mpa_version = db_version,
        mpa_taxProfile = mtprofile
      )
    }
  )
  return(result)
}


#' Join multiple MetaPhlAn profiles
#'
#' Combines multiple `metaphlanProfile` objects into a single profile by joining
#' their taxonomic abundance tables on the `clade_name` column.
#'
#' @param ... Multiple objects of class `metaphlanProfile` to merge.
#'
#' @return An object of class `metaphlanProfile` containing the merged taxonomic profiles.
#'   All input profiles must use the same MetaPhlAn database version.
#'
#' @details
#' - Non-`metaphlanProfile` arguments are silently dropped with a warning.
#' - Profiles from different MetaPhlAn database versions will cause an error.
#' - Numeric suffixes will automatically be added to duplicate samples (e.g. Sample_1 -> Sample_1.1, Sample_1.2).
#' - Missing values (NA) after joining are replaced with 0.0.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mtp1 <- load_metaphlan("path/to/profile1.txt")
#' mtp2 <- load_metaphlan("path/to/profile2.txt")
#' merged_profile <- join_mpa_profiles(mtp1, mtp2)
#' }
join_mpa_profiles <- function(...) {
  arglist <- list(...)
  mpa_profiles <- Filter(function(x) methods::is(x, "metaphlanProfile"), arglist)

  if (length(mpa_profiles) != length(arglist)) {
    warning('At least one input object is not of class "metaphlanProfile" and was dropped.')
  }

  db_versions <- lapply(mpa_profiles, function(x) x@mpa_version)

  if (!all(vapply(db_versions, length, integer(1L)) > 0)) {
    stop("One or more profiles have an empty mpa_version (character(0)).")
  }

  all_same_db <- length(unique(unlist(db_versions))) == 1
  # all_same_db <- length(unique(unlist(lapply(mpa_profiles, function(x) x@mpa_version)))) == 1
  if (!all_same_db) {
    stop(
      "Profiles from different versions of MetaPhlAn were found!\n",
      "Please make sure the used 'mpa' database is always the same."
    )
  }

  db_version <- mpa_profiles[[1]]@mpa_version

  mpa_profiles <- lapply(mpa_profiles, function(x) x@mpa_taxProfile)

  # Get sample columns (exclude clade_name)
  sample_cols <- lapply(mpa_profiles, function(df) setdiff(names(df), "clade_name"))

  # Find overlapping sample names
  all_samples <- unlist(sample_cols)
  duplicates <- names(which(table(all_samples) > 1))

  if (length(duplicates) > 0) {
    warning(
      "Duplicate sample names detected: ", paste(duplicates, collapse = ", "),
      ". Automatically disambiguating with suffixes."
    )
  }

  # Disambiguate column names by adding suffixes (e.g. Sample_1 -> Sample_1.1, Sample_1.2)
  mpa_profiles <- mapply(function(df, suffix) {
    conflict_cols <- intersect(duplicates, names(df))
    if (length(conflict_cols) > 0) {
      colnames(df)[colnames(df) %in% conflict_cols] <- paste0(colnames(df)[colnames(df) %in% conflict_cols], ".", suffix)
    }
    df
  }, mpa_profiles, seq_along(mpa_profiles), SIMPLIFY = FALSE)

  # Join profiles
  mpa_profiles_df <- tryCatch(
    {
      purrr::reduce(mpa_profiles, dplyr::full_join, by = "clade_name")
    },
    error = function(condition) {
      stop('Please ensure that the MetaPhlAn profiles have a "clade_name" column to join on.')
    }
  )

  mpa_profiles_df[is.na(mpa_profiles_df)] <- 0.0

  result <- new("metaphlanProfile",
    mpa_version = db_version,
    mpa_taxProfile = mpa_profiles_df
  )
  return(result)
}
