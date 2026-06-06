#' A S4 class extending phyloseq for microbiome data analysis
#'
#' This class inherits from \code{\link[phyloseq]{phyloseq}} and adds a slot for
#' additional information, such as lowest taxonomic rank, MetaPhlAn database
#' version, applied transformations and/or filters.
#'
#' @slot info A `list` containing metadata.
#' @slot otu_table A `matrix` of class `otu_table` with taxa as rows and samples as columns.
#' @slot tax_table A named `character matrix` of class `taxonomyTable`, with taxon IDs as row names and taxonomic ranks as columns.
#' @slot sam_data An object of class `sample_data` containing sample metadata.
#' @slot phy_tree An object of class `phylo` representing a phylogenetic tree.
#' @slot refseq An object of a class inheriting from \code{XStringSet-class} (e.g., \code{DNAStringSet}) containing biological sequences. See \code{\link[Biostrings]{XStringSet}}.
#' @importClassesFrom phyloseq phyloseq
#' @export
setClass("microEDA",
  contains = "phyloseq",
  slots = list(info = "list"),
)


#### microEDA replacement methods for phyloseq slots ####
# Implemented to preserve the microEDA class instead of coercing to phyloseq
### //////////////////////////////////////////////////////////////////////////###
#' Replace OTU table in a microEDA object
#'
#' @param x A `microEDA` object.
#' @param value An `otu_table` object (or coercible to one).
#' @return The updated `microEDA` object.
#' @method otu_table<- microEDA
#' @aliases otu_table<-,microEDA,otu_table-method
#' @importFrom phyloseq otu_table<-
#' @export
setMethod(
  "otu_table<-", c("microEDA", "otu_table"),
  function(x, value) {
    if (!is(value, "otu_table")) {
      value <- as(value, "otu_table")
    }
    x@otu_table <- value
    validObject(x)
    x
  }
)


#' Replace taxonomic table in a microEDA object
#'
#' @param x A `microEDA` object.
#' @param value A `taxonomyTable` object (or coercible to one).
#' @return The updated `microEDA` object.
#' @method tax_table<- microEDA
#' @aliases tax_table<-,microEDA,taxonomyTable-method
#' @importFrom phyloseq tax_table<-
#' @export
setMethod(
  "tax_table<-", c("microEDA", "taxonomyTable"),
  function(x, value) {
    if (!is(value, "taxonomyTable")) {
      value <- as(value, "taxonomyTable")
    }
    x@tax_table <- value
    validObject(x)
    x
  }
)


#' Replace sample data in a microEDA object
#'
#' @param x A `microEDA` object.
#' @param value A `sample_data` object (or coercible to one).
#' @return The updated `microEDA` object.
#' @method sample_data<- microEDA
#' @aliases sample_data<-,microEDA,sample_data-method
#' @importFrom phyloseq sample_data
setGeneric("sample_data<-", function(x, value) standardGeneric("sample_data<-"))
#' @export
setMethod(
  "sample_data<-", c("microEDA", "sample_data"),
  function(x, value) {
    if (!is(value, "sample_data")) {
      value <- as(value, "sample_data")
    }
    x@sam_data <- value
    validObject(x)
    x
  }
)


#' Replace phylogenetic tree in a microEDA object
#'
#' @param x A `microEDA` object.
#' @param value A `phylo` object.
#' @return The updated `microEDA` object.
#' @method phy_tree<- microEDA
#' @aliases phy_tree<-,microEDA,phylo-method
#' @importFrom phyloseq phy_tree<-
#' @export
setMethod(
  "phy_tree<-", c("microEDA", "phylo"),
  function(x, value) {
    if (!is(value, "phylo")) {
      stop("value must be a phylo object")
    }
    x@phy_tree <- value
    validObject(x)
    x
  }
)


#### microEDA methods for info slot ####
### //////////////////////////////////////////////////////////////////////////###
.valid_keys <- c("taxrank", "transforms", "filter_history", "mpa_version", "filtered_taxa")

#' Get or set the info slot from a microEDA object
#'
#' @description
#' Accessor methods for either retrieving or setting the info slot in a
#' \linkS4class{microEDA} object.
#' Additionally, each expected element of the info slot can be retrieved or set
#' individually.
#' Valid keys for the `microEDA` info slot can be displayed with `infoKeys()`,
#' while filled info fields of a `microEDA` object can be displayed with
#' `infoFields()`.
#'
#' @param object A \linkS4class{microEDA} object.
#' @param value Only used in replacement methods (setters).
#' \describe{
#'  \item{For `info<-`: }{A `list` containing metadata. Valid metadata keys are: 'taxrank', 'transforms', 'filter_history', 'mpa_version'.}
#'  \item{For `taxrank<-`: }{A `character` string specifying the taxonomic rank (e.g., "Species", "Genus").}
#'  \item{For `transforms<-`: }{A `character` vector of applied transformation(s) (e.g., c("TSS", "log2")).}
#'  \item{For `mpa_version<-`: }{A `character` string specifying the MetaPhlAn version (e.g., "#mpa_vJun23_CHOCOPhlAnSGB_202307").}
#'  \item{For `filtered_taxa<-`: }{A `list` with otu_table and tax_table (as matrices) of features that were filtered out.}
#'  \item{For `filter_history<-`: }{A named list specifying used filter parameters. Expected elements include:
#'    \itemize{
#'      \item \code{min_abundance}: numeric threshold for minimum abundance.
#'      \item \code{min_prevalence}: numeric threshold for minimum prevalence.
#'      \item \code{group_var}: optional grouping variable for stratified filtering.
#'      \item \code{abundance_criterion}: filtering criterion applied to abundance values.
#'      \item \code{group_requirement}: filter passing requirement for group-wise filtering.
#'      \item \code{keep_filtered}: logical indicating whether filtered-out features were preserved in `@info$filtered_taxa`.
#'      \item \code{n_removed}: integer indicating the number of features removed.
#'      }
#'    See \code{\link{filter_features}} for details on filtering arguments.
#'    }
#'  }
#' @return \describe{
#'  \item{For `info()`: }{A `list` containing metadata.}
#'  \item{For `infoKeys()`: }{Lists valid keys for `microEDA` info slot.}
#'  \item{For `infoFields()`: }{Lists available field(s) in `object` info slot.}
#'  \item{For `taxrank()`: }{A `character` string with the lowest taxonomic rank.}
#'  \item{For `transforms()`: }{A `character` vector of applied transformations.}
#'  \item{For `mpa_version()`: }{A `character` string with the used the MetaPhlAn database version.}
#'  \item{For `filtered_taxa()`: }{A `list` with otu_table and tax_table of features that were filtered out.}
#'  \item{For `filter_history()`: }{A `list` of used filter parameters.}
#' }
#' Replacement methods (`<-`) will return a `microEDA` object with the updated
#' corresponding slots.
#' @details
#' `taxrank`: represents the lowest available taxonomic rank.
#'
#' `transforms`: represents transformation(s) that were applied to the abundance data (otu_table).
#'
#' `mpa_version`: represents the MetaPhlAn database version that was used to create the profile.
#'
#' `filtered_taxa`: represents the otu_table and tax_table of features that were filtered out.
#'
#' `filter_history`: represents filters that were applied to the abundance data (otu_table).
#'
#' @rdname info-accessors
#' @aliases info,microEDA-method
#' @docType methods
#' @export
setGeneric("info", function(object) standardGeneric("info"))

setMethod(
  "info", "microEDA",
  function(object) {
    return(object@info)
  }
)


#' @rdname info-accessors
#' @aliases info<-,microEDA-method
#' @export
setGeneric("info<-", function(object, value) standardGeneric("info<-"))

setMethod("info<-", "microEDA", function(object, value) {
  if (!is.list(value)) stop("'info' must be a named list.")

  if (any(!names(value) %in% .valid_keys)) {
    invalid <- names(value)[!names(value) %in% .valid_keys]
    warning(
      "Unknown metadata added to 'info': ", paste(invalid, collapse = ", "),
      ". Expected keys are: ", paste(.valid_keys, collapse = ", ")
    )
  }
  object@info <- value
  return(object)
})


#' @rdname info-accessors
#' @export
infoKeys <- function() {
  message(
    "Expected keys for 'microEDA' info slot are: ", paste(.valid_keys, collapse = ", ")
  )
}

#' @rdname info-accessors
#' @aliases infoFields,microEDA-method
#' @export
setGeneric("infoFields", function(object) standardGeneric("infoFields"))
setMethod("infoFields", "microEDA", function(object) {
  fields <- names(object@info)[!vapply(object@info, is.null, logical(1L))]
  message(
    "Available field(s) in '", deparse(substitute(object)), "' info slot are: ",
    paste(fields, collapse = ", ")
  )
})


#' @rdname info-accessors
#' @aliases taxrank,microEDA-method
#' @export
setGeneric("taxrank", function(object) standardGeneric("taxrank"))
setMethod("taxrank", "microEDA", function(object) object@info[["taxrank"]])

#' @rdname info-accessors
#' @aliases taxrank<-,microEDA-method
#' @export
setGeneric("taxrank<-", function(object, value) standardGeneric("taxrank<-"))
setMethod("taxrank<-", "microEDA", function(object, value) {
  object@info[["taxrank"]] <- value
  validObject(object)
  object
})


#' @rdname info-accessors
#' @aliases transforms,microEDA-method
#' @export
setGeneric("transforms", function(object) standardGeneric("transforms"))
setMethod("transforms", "microEDA", function(object) object@info[["transforms"]])

#' @rdname info-accessors
#' @aliases transforms<-,microEDA-method
#' @export
setGeneric("transforms<-", function(object, value) standardGeneric("transforms<-"))
setMethod("transforms<-", "microEDA", function(object, value) {
  object@info[["transforms"]] <- value
  validObject(object)
  object
})


#' @rdname info-accessors
#' @aliases mpa_version,microEDA-method
#' @export
#' @seealso \code{\link{mpa_version,metaphlanProfile-method}}
setMethod("mpa_version", "microEDA", function(object) object@info[["mpa_version"]])

#' @rdname info-accessors
#' @aliases mpa_version<-,microEDA-method
#' @export
setMethod("mpa_version<-", "microEDA", function(object, value) {
  object@info[["mpa_version"]] <- value
  validObject(object)
  object
})


#' @rdname info-accessors
#' @aliases filtered_taxa,microEDA-method
#' @export
setGeneric("filtered_taxa", function(object) standardGeneric("filtered_taxa"))
setMethod("filtered_taxa", "microEDA", function(object) object@info[["filtered_taxa"]])

#' @rdname info-accessors
#' @aliases filtered_taxa<-,microEDA-method
#' @export
setGeneric("filtered_taxa<-", function(object, value) standardGeneric("filtered_taxa<-"))
setMethod("filtered_taxa<-", "microEDA", function(object, value) {
  # if (!inherits(value, "phyloseq")) stop("'value' must be a phyloseq object containing otu and tax table of filtered out taxa.")
  if (!is.list(value)) stop("'value' must be a named list containing otu and tax table of filtered out taxa.")
  object@info[["filtered_taxa"]] <- value
  validObject(object)
  object
})


#' @rdname info-accessors
#' @aliases filter_history,microEDA-method
#' @export
setGeneric("filter_history", function(object) standardGeneric("filter_history"))
setMethod("filter_history", "microEDA", function(object) object@info[["filter_history"]])

#' @rdname info-accessors
#' @aliases filter_history<-,microEDA-method
#' @export
setGeneric("filter_history<-", function(object, value) standardGeneric("filter_history<-"))
setMethod("filter_history<-", "microEDA", function(object, value) {
  .valid_filter_arguments <- c(
    "min_abundance",
    "min_prevalence",
    "group_var",
    "abundance_criterion",
    "group_requirement",
    "keep_filtered",
    "n_removed"
  )

  if (!is.list(value)) stop("'filter_history' must be a named list.")

  if (any(!names(value) %in% .valid_filter_arguments)) {
    invalid <- names(value)[!names(value) %in% .valid_filter_arguments]
    stop(
      "Unknown filter arguments: ", paste(invalid, collapse = ", "),
      ". Expected keys are: ", paste(.valid_keys, collapse = ", ")
    )
  }

  if (is.null(value$group_var)) {
    value$group_var <- ""
    value$group_requirement <- ""
  }

  # Initialize filter_history if not present
  if (is.null(object@info[["filter_history"]])) {
    object@info[["filter_history"]] <- value
  } else { # Append current filter run
    object@info[["filter_history"]] <- mapply(c, object@info[["filter_history"]], value, SIMPLIFY = FALSE)
  }

  validObject(object)
  object
})


#### Methods for MetaPhlAn profiles ####
### //////////////////////////////////////////////////////////////////////// ###

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

  tax_tab <- .mpa_prepare_tax_tab(tax_profile@mpa_taxProfile)
  abund_tab <- .mpa_prepare_abundance_tab(tax_profile@mpa_taxProfile, tax_tab)
  # Drop no longer needed clade_name column
  tax_tab <- tax_tab[, !colnames(tax_tab) %in% c("clade_name")]

  if (any(abund_tab < 0, na.rm = TRUE)) .show_error("OTU table contains negative values.")

  if (.is_proportion(abund_tab, silent = TRUE)) {
    transforms <- "TSS"
  } else if (.is_counts(abund_tab, silent = TRUE)) {
    transforms <- NULL
  } else {
    stop("'tax_profile' contains neither counts nor relative abundances - data was likely transformed.")
  }

  if (is.null(metadata)) {
    microbiome_exp <- new(
      "microEDA",
      otu_table = phyloseq::otu_table(abund_tab, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_tab),
      sam_data = NULL,
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_tab), n = 1),
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
      otu_table = phyloseq::otu_table(abund_tab, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_tab),
      sam_data = phyloseq::sample_data(metadata),
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_tab), n = 1),
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
      otu_table = phyloseq::otu_table(abund_tab, taxa_are_rows = TRUE),
      tax_table = phyloseq::tax_table(tax_tab),
      sam_data = NULL,
      phy_tree = NULL,
      refseq = NULL,
      info = list(
        taxrank = tail(colnames(tax_tab), n = 1),
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


#### Methods for phyloseq objects ####
### //////////////////////////////////////////////////////////////////////// ###

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

  abund_tab <- phyloseq::otu_table(tax_profile)
  tax_tab <- tax_profile@tax_table

  # Ensure abundance table is in expected orientation
  if (!phyloseq::taxa_are_rows(abund_tab)) {
    abund_tab <- phyloseq::t(abund_tab)
  }

  if (any(abund_tab < 0, na.rm = TRUE)) .show_error("OTU table contains negative values.")

  if (is.null(tax_tab)) .show_error('The "tax_table" slot must not be empty for "phyloseq" objects.')

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

  if (.is_proportion(abund_tab, silent = TRUE)) {
    transforms <- "TSS"
  } else if (.is_counts(abund_tab, silent = TRUE)) {
    transforms <- NULL
  } else {
    stop("otu_table of 'tax_profile' contains neither counts nor relative abundances - data was likely transformed.")
  }

  microbiome_exp <- new(
    "microEDA",
    otu_table = abund_tab,
    tax_table = tax_tab,
    sam_data = metadata,
    phy_tree = tax_profile@phy_tree,
    refseq = tax_profile@refseq,
    info = list(
      taxrank = tail(colnames(tax_tab), n = 1),
      transforms = transforms
    )
  )

  return(microbiome_exp)
}

### //////////////////////////////////////////////////////////////////////// ###
#### Methods for microEDA objects ####

#' Constructor for microEDA objects
#'
#' @param tax_profile An object of class `metaphlanProfile` or `phyloseq`.
#' @param metadata Optional `data.frame` with sample metadata.
#' @param sample_column `Character` string specifying the column name in
#'  `metadata` that contains sample identifiers (e.g., "Samples").
#' @return An object of class `microEDA`.
#' @rdname microEDA-class
#' @docType methods
#' @exportMethod microEDA
#' @details
#' When passed a \linkS4class{metaphlanProfile}, the profile will be returned at
#' the lowest available taxonomic rank, which usually is at strain rank. To
#' agglomerate the profile at a higher rank, please see [`agglomerate_taxa()`].
#'
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


#' Show method for microEDA objects
#'
#' Displays a summary of a `microEDA` object, including the standard
#' `phyloseq` output and custom information from the `info` slot.
#'
#' @param object A `microEDA` object.
#' @exportMethod show
#' @rdname microEDA-class
setMethod("show", "microEDA", function(object) {
  cat("=== microEDA - an extended phyloseq object ===\n\n")

  # Call the parent phyloseq show method
  callNextMethod()

  # Print custom 'info' slot header
  cat("\n--- microEDA Info ---\n")

  # Display info slot contents
  if (length(object@info) == 0) {
    cat("info: <empty>\n")
  } else {
    for (field_name in names(object@info)) {
      val <- object@info[[field_name]]

      if (is.null(val)) {
        cat(sprintf("%s: <NULL>\n", field_name))
      } else if (is.vector(val) && length(val) == 1) {
        cat(sprintf("%s: %s\n", field_name, val))
      } else if (field_name == "filter_history" && is.list(val) && length(val) > 0) {
        # Special handling for filter_history: count the depth of steps
        # Assumes all steps have the same length/structure. Modify if
        # corresponding section in filter_features changes
        n_steps <- length(val[[1]])
        step_word <- if (n_steps == 1) "step" else "steps"
        cat(sprintf("%s: [ %d filter %s applied ]\n", field_name, n_steps, step_word))
      } else if (is.list(val) && length(val) > 5) {
        # General fallback for other long lists
        cat(sprintf("%s: <list of %d items>\n", field_name, length(val)))
      } else if (field_name == "filtered_taxa" && is.list(val)) {
        # Special handling for filtered_taxa
        otu_mat <- val$filtered_otu
        tax_mat <- val$filtered_tax

        # Ensure they are matrices before accessing dimensions
        if (is.matrix(otu_mat) && is.matrix(tax_mat)) {
          n_taxa <- nrow(otu_mat)
          n_samples <- ncol(otu_mat)
          n_ranks <- ncol(tax_mat)

          cat(sprintf(
            "%s: [ %d taxa by %d taxonomic ranks in %d samples ]\n",
            field_name, n_taxa, n_ranks, n_samples
          ))
        } else {
          # Fallback if structure is unexpected
          cat(sprintf("%s: <invalid structure>\n", field_name))
        }
      } else if (inherits(val, "phyloseq")) {
        cat(sprintf(
          "%s: phyloseq-object [%d taxa, %d samples]\n",
          field_name, phyloseq::ntaxa(val), phyloseq::nsamples(val)
        ))
      } else {
        cat(sprintf("%s: %s\n", field_name, paste(val, collapse = ", ")))
      }
    }
  }
})
