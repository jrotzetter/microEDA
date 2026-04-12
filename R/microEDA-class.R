#' An S4 class extending phyloseq for microbiome data analysis
#'
#' This class inherits from phyloseq and adds a slot for additional information,
#' such as lowest taxonomic rank, MetaPhlAn database version, applied
#' transformations and/or filters.
#'
#' @slot info A list containing metadata.
#' @importClassesFrom phyloseq phyloseq
#' @export
setClass("microEDA",
  contains = "phyloseq",
  slots = list(info = "list"),
)


#' Extract taxonomic table
#'
#' @param object An object of class \code{microEDA}.
#' @return The tax_table component of the microEDA/phyloseq object.
#' @rdname getTaxTable-methods
#' @docType methods
#' @exportMethod getTaxTable
setGeneric("getTaxTable", function(object) standardGeneric("getTaxTable"))
#' @rdname getTaxTable-methods
#' @aliases getTaxTable,microEDA-method
setMethod("getTaxTable", "microEDA", function(object) {
  object@tax_table
})


#' Extract sample data
#'
#' @param object An object of class \code{microEDA}.
#' @return The sam_data component of the microEDA/phyloseq object.
#' @rdname getSampleData-methods
#' @docType methods
#' @exportMethod getSampleData
setGeneric("getSampleData", function(object) standardGeneric("getSampleData"))
#' @rdname getSampleData-methods
#' @aliases getSampleData,microEDA-method
setMethod("getSampleData", "microEDA", function(object) {
  object@sam_data
})


#' Check taxonomic rank of the profile
#'
#' @param object An object of class \code{microEDA}.
#' @return The taxonomic rank stored in the \code{info} slot.
#' @rdname checkProfileRank-methods
#' @docType methods
#' @exportMethod checkProfileRank
setGeneric("checkProfileRank", function(object) standardGeneric("checkProfileRank"))
#' @rdname checkProfileRank-methods
#' @aliases checkProfileRank,microEDA-method
setMethod("checkProfileRank", "microEDA", function(object) {
  object@info$taxrank
})


### //////////////////////////////////////////////////////////////////////// ###
#### Methods for MetaPhlAn profiles ####

#' Constructor for microEDA objects from MetaPhlAn profiles
#'
#' Internal function to create a `microEDA` object by combining taxonomic
#' abundance data and optional sample metadata.
#'
#' @param tax_profile A `metaphlanProfile` object.
#' @param metadata Optional `data.frame` containing sample metadata.
#' @param sample_column `Character` string specifying the column name in
#'  `metadata` that contains sample identifiers (e.g., "Samples").
#' @importFrom utils tail
#'
#' @return An object of class `microEDA`.
#' @keywords internal
#' @noRd
.metaphlanConstructor <- function(tax_profile, metadata = NULL, sample_column = NULL) {
  stopifnot(inherits(tax_profile, "metaphlanProfile"))

  tax_table <- .mpa_prepare_tax_tab(tax_profile@mpa_taxProfile)
  abund_table <- .mpa_prepare_abundance_tab(tax_profile@mpa_taxProfile, tax_table)
  # Drop no longer needed clade_name column
  tax_table <- tax_table[, !colnames(tax_table) %in% c("clade_name")]

  if (any(abund_table < 0, na.rm = TRUE)) .show_error("OTU table contains negative values.")

  if (!.is_proportion(abund_table)) {
    transforms <- NULL
  } else {
    transforms <- "relabund"
  }

  if (is.null(metadata)) {
    microbiome_exp <- new(
      "microEDA",
      otu_table = phyloseq::otu_table(abund_table, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_table),
      sam_data = NULL,
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_table), n = 1),
        transforms = transforms,
        mpa_version = tax_profile@mpa_version
      )
    )
  } else if (.mpa_check_names(tax_profile@mpa_taxProfile, metadata, sample_column)) {
    rownames(metadata) <- NULL
    metadata <- tibble::column_to_rownames(metadata, sample_column)

    # Make sure no reserved column names are in the metadata
    metadata <- .check_var_names(metadata)

    microbiome_exp <- new(
      "microEDA",
      otu_table = phyloseq::otu_table(abund_table, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_table),
      sam_data = phyloseq::sample_data(metadata),
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_table), n = 1),
        transforms = transforms,
        mpa_version = tax_profile@mpa_version
      )
    )
  } else {
    .show_warning(paste(
      'Sample names between "tax_profile" and "metadata"',
      'did not match. Ignoring "metadata".'
    ))

    microbiome_exp <- new(
      "microEDA",
      otu_table = phyloseq::otu_table(abund_table, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_table),
      sam_data = NULL,
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_table), n = 1),
        transforms = transforms,
        mpa_version = tax_profile@mpa_version
      )
    )
  }
  return(microbiome_exp)
}


#' Check if sample names in MetaPhlAn profile and metadata match
#'
#' @param mtphlan_profile `mpa_taxProfile` from a metaphlanProfile object.
#' @param metadata `data.frame` with sample metadata.
#' @param sample_col column in metadata to use as sample names.
#' @return `logical` indicating whether all sample names match.
#' @keywords internal
#' @noRd
.mpa_check_names <- function(mtphlan_profile, metadata, sample_col) {
  stopifnot(is.data.frame(mtphlan_profile))
  # Remove clade_name and all other non-numeric columns
  mtphlan_profile <- mtphlan_profile |> dplyr::select(dplyr::where(is.numeric))

  mp_samples <- names(mtphlan_profile)

  mt_samples <- metadata[[sample_col]]

  names_ok <- TRUE

  if (!all(mp_samples %in% mt_samples)) {
    mp_missing_index <- which(!mp_samples %in% mt_samples)
    mp_missing <- names(mtphlan_profile[mp_missing_index])
    warning("Samples ", paste0('"', mp_missing, '"', collapse = ", "), " of the MetaPhlAn profile are not in the metadata!")
    names_ok <- FALSE
  }

  if (!all(mt_samples %in% mp_samples)) {
    mt_missing_index <- which(!mt_samples %in% mp_samples)
    mt_missing <- metadata[mt_missing_index, sample_col]
    warning("Metadata samples ", paste0('"', mp_missing, '"', collapse = ", "), " are not in the MetaPhlAn profile!")
    names_ok <- FALSE
  }
  return(names_ok)
}


#' Prepare taxonomic table from MetaPhlAn profile
#'
#' Internal function to convert a MetaPhlAn profile into a standardized
#' taxonomic table with proper rank prefixes and structure.
#'
#' @param metaphlan_profile A data frame containing a `clade_name` column
#'   with full lineage strings (e.g., "k__Bacteria|p__Firmicutes|...").
#' @return A matrix or data frame where rows are taxa, columns are taxonomic
#'   ranks (e.g., "Phylum", "Class"), and the first column is `clade_name`.
#'   Rows with unclassified lowest ranks are handled appropriately.
#' @keywords internal
#' @noRd
.mpa_prepare_tax_tab <- function(metaphlan_profile) {
  metaphlan_profile <- metaphlan_profile$clade_name
  split_clades <- strsplit(metaphlan_profile, "|", fixed = TRUE)
  max_ranks <- max(lengths(split_clades))

  # Get all lowest taxa, then retrieve the prefix
  last_taxa <- vapply(split_clades, tail, n = 1, character(1L))
  tax_chars <- unique(substr(last_taxa, 1, 1))

  if (!length(tax_chars) > 1) {
    stop("Only one taxonomic rank found in MetaPhlAn profile. At least two ranks are required.")
  }

  # Retrieve the taxonomic levels present in the data
  tax_cols <- names(.tax_ranks)[!is.na(match(.tax_ranks, tax_chars))]

  # Set length of each vector to the maximum length, automatically filling
  # missing positions with NA
  padded_data <- lapply(split_clades, `length<-`, max_ranks)
  # Remove any rank prefixes
  padded_data <- lapply(padded_data, .trim_rank_prefix)
  taxtab <- do.call(rbind, padded_data)

  # Add taxonomic levels as column names
  colnames(taxtab) <- tax_cols[1:ncol(taxtab)]

  # Check if last taxon starts either with the prefix of the lowest rank or
  # UNCLASSIFIED...
  pattern <- paste0("^", .tax_ranks[max_ranks], "__|^UNCLASSIFIED")
  lowest_tax_ind <- grepl(pattern, last_taxa)
  # ...then subset the tax table and metaphlan profile to only include these rows
  taxtab <- taxtab[lowest_tax_ind, ]
  metaphlan_profile <- metaphlan_profile[lowest_tax_ind]

  # Insert clade_name at lowest taxonomic level at beginning
  taxtab <- cbind(
    clade_name = metaphlan_profile,
    taxtab
  )

  if (any(taxtab[, "clade_name"] == "UNCLASSIFIED")) {
    taxtab[taxtab[, "clade_name"] == "UNCLASSIFIED", ] <- "UNCLASSIFIED"
  }
  rownames(taxtab) <- taxtab[, ncol(taxtab)]
  return(taxtab)
}


#' Prepare abundance table from MetaPhlAn profile
#'
#' Internal function to subset and format the abundance data matrix from a
#' MetaPhlAn profile, aligning it with the taxonomic table.
#'
#' @param metaphlan_profile A data frame containing abundance values and a
#'   `clade_name` column.
#' @param tax_tab A taxonomic table (e.g., from `.mpa_prepare_tax_tab`) with
#'   row names corresponding to clade names of the lowest taxonomic rank.
#' @return A matrix of abundance values with row names matching `tax_tab`.
#' @keywords internal
#' @noRd
.mpa_prepare_abundance_tab <- function(metaphlan_profile, tax_tab) {
  # Subset metaphlan_profile to only include rows of the lowest taxonomic rank,
  # before dropping the clade_name column
  metaphlan_profile <- metaphlan_profile[metaphlan_profile$clade_name %in% tax_tab[, "clade_name"], ]
  metaphlan_profile <- metaphlan_profile[, !names(metaphlan_profile) %in% c("clade_name")]
  rownames(metaphlan_profile) <- rownames(tax_tab)
  metaphlan_profile <- as(metaphlan_profile, "matrix")
  return(metaphlan_profile)
}


### //////////////////////////////////////////////////////////////////////// ###
#### Methods for phyloseq objects ####

#' Constructor for microEDA objects from a phyloseq object
#'
#' Internal function to create a microEDA object by combining an otu_table
#' and optional sample metadata (sam_data).
#'
#' @param tax_profile A `phyloseq` object.
#' @param metadata Optional `data.frame` containing sample metadata.
#' @param sample_column `Character` string specifying the column name in
#'  `metadata` that contains sample identifiers (e.g., "Samples").
#' @importFrom utils tail
#'
#' @return An object of class `microEDA`.
#' @keywords internal
#' @noRd
.phyloseqConstructor <- function(tax_profile, metadata = NULL, sample_column = NULL) {
  stopifnot(inherits(tax_profile, "phyloseq"))

  abund_table <- phyloseq::otu_table(tax_profile)
  tax_table <- tax_profile@tax_table

  # Ensure abundance table is in expected orientation
  if (!phyloseq::taxa_are_rows(abund_table)) {
    abund_table <- phyloseq::t(abund_table)
  }

  if (any(abund_table < 0, na.rm = TRUE)) .show_error("OTU table contains negative values.")

  if (is.null(tax_table)) .show_error('The "tax_table" slot must not be empty for "phyloseq" objects.')

  if (!is.null(sample_column) && is.null(metadata)) .show_error('"metadata" must be provided when "sample_column" is specified to add metadata.')

  # If metadata is provided, use sample_column to set row names
  if (!is.null(metadata)) {
    if (is.null(sample_column)) {
      .show_error('"sample_column" must be specified when providing metadata.')
    }

    if (!(sample_column %in% colnames(metadata))) {
      .show_error(paste0('Column "', sample_column, '" not found in metadata.'))
    }

    # Extract sample names from phyloseq object for checks and ordering
    ps_sample_names <- phyloseq::sample_names(tax_profile)

    # Set sample names from specified column
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- NULL
    metadata <- tibble::column_to_rownames(metadata, sample_column)
    metadata <- phyloseq::sample_data(metadata)

    # Subset and reorder metadata to match phyloseq sample order
    if (!all(ps_sample_names %in% rownames(metadata))) {
      .show_error("Some samples in phyloseq object are missing from metadata.")
    }
    # Reorder and ensure all samples present
    metadata <- metadata[ps_sample_names, , drop = FALSE]
  } else {
    # If no metadata provided, extract from phyloseq
    metadata <- phyloseq::sample_data(tax_profile, errorIfNULL = FALSE)
  }

  # Make sure no reserved column names are in the metadata
  metadata <- .check_var_names(metadata)

  if (!.is_proportion(abund_table, silent = TRUE)) {
    transforms <- NULL
  } else {
    transforms <- "relabund"
  }

  microbiome_exp <- new(
    "microEDA",
    otu_table = abund_table,
    tax_table = tax_table,
    sam_data = metadata,
    phy_tree = tax_profile@phy_tree,
    refseq = tax_profile@refseq,
    info = list(
      taxrank = tail(colnames(tax_table), n = 1),
      transforms = transforms
    )
  )

  return(microbiome_exp)
}

### //////////////////////////////////////////////////////////////////////// ###
#### Methods for microEDA objects ####

#' Constructor for microEDA objects
#'
#' @param tax_profile An object of class "metaphlanProfile" or "phyloseq".
#' @param metadata Optional data.frame with sample metadata.
#' @param sample_column `Character` string specifying the column name in
#'  `metadata` that contains sample identifiers (e.g., "Samples").
#' @return An object of class `microEDA`.
#' @rdname microEDA-class
#' @docType methods
#' @exportMethod microEDA
#'
#' @examples
#' data("GlobalPatterns", package = "phyloseq")
#' me <- microEDA(GlobalPatterns)
setGeneric("microEDA", function(tax_profile,
                                metadata = NULL,
                                sample_column = NULL) {
  if (!inherits(tax_profile, "phyloseq") && !inherits(tax_profile, "metaphlanProfile")) {
    stop("'tax_profile' must be a metaphlanProfile or phyloseq object", call. = FALSE)
  }
  standardGeneric("microEDA")
})


#' @rdname microEDA-class
#' @aliases microEDA,metaphlanProfile-method
setMethod("microEDA", "metaphlanProfile", function(tax_profile,
                                                   metadata = NULL,
                                                   sample_column = NULL) {
  .metaphlanConstructor(tax_profile,
    metadata = metadata,
    sample_column = sample_column
  )
})


#' @rdname microEDA-class
#' @aliases microEDA,phyloseq-method
setMethod("microEDA", "phyloseq", function(tax_profile,
                                           metadata = NULL,
                                           sample_column = NULL) {
  .phyloseqConstructor(tax_profile,
    metadata = metadata,
    sample_column = sample_column
  )
})
