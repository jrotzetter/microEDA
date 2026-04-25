#' Fill unclassified taxonomic rank using last known parent clade
#'
#' This internal function replaces missing or empty taxonomic entries (including `NA`, `""`, `" "`, `"\t"`)
#' with "Unclassified". For the specified target rank, it further relabels "Unclassified" entries
#' as "Unclassified {Last_Known_Parent}" using the last known parent clade.
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#'   (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species) in hierarchical order,
#'   and rows represent taxa.
#' @param taxrank A `character` string specifying the target taxonomic rank (e.g., "Genus")
#'   to process and relabel.
#' @param silent A `logical` value indicating whether to suppress warnings.
#'   If `FALSE`, a warning is issued if any rows are fully unclassified across all ranks.
#'
#' @return A `matrix` with:
#'   - Missing/empty values replaced by "Unclassified".
#'   - Only columns from the first up to and including `taxrank`.
#'   - "Unclassified" entries in the `taxrank` column relabeled as "Unclassified {Parent}"
#'     where `{Parent}` is the last non-"Unclassified" value in the row.
#'   - Fully unclassified rows remain as "Unclassified".
#'
#' @examples
#' taxa <- data.frame(
#'   Kingdom = c("Kingdom1", "Kingdom2", "Kingdom3", "Kingdom4", "Kingdom5", "Kingdom6", ""),
#'   Phylum = c("Phylum1", "", "Phylum3", "Phylum4", "Phylum5", "", NA),
#'   Class = c("Class1", "Class2", " ", "\t", NA, "Class6", ""),
#'   Order = c("Order1", "Order2", NA, "Order4", NA, "", " ")
#' )
#' rownames(taxa) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5", "OTU6", "OTU7")
#'
#' .fill_unclassified_rank(taxa, "Order")
#'
#' @keywords internal
#' @noRd
.fill_unclassified_rank <- function(taxa, taxrank, silent = TRUE) {
  stopifnot(is.data.frame(taxa) || is.matrix(taxa))
  stopifnot(is.character(taxrank))
  stopifnot(taxrank %in% colnames(taxa))

  # Convert to character matrix for speed
  taxa <- as.matrix(taxa)

  # Replace all NA's or empty strings
  empty_vals <- c(NA, "", " ", "\t")
  taxa[taxa %in% empty_vals] <- "Unclassified"

  # Keep only columns up to and including 'taxrank'
  idx <- which(colnames(taxa) == taxrank)
  taxa <- taxa[, 1:idx, drop = FALSE]

  # Find rows where target rank is "Unclassified"
  is_unclass <- taxa[, idx] == "Unclassified"
  if (!any(is_unclass)) {
    return(taxa)
  }

  # Find last non-"Unclassified" value for each row
  parent_ranks <- taxa[is_unclass, 1:(idx - 1), drop = FALSE]
  last_valid <- vapply(
    seq_len(nrow(parent_ranks)),
    function(i) {
      row <- parent_ranks[i, , drop = FALSE]
      valid <- tail(row[row != "Unclassified"], 1)
      if (length(valid) == 0) "" else valid
    },
    character(1L)
  )

  if (!silent) {
    n_all_unclass <- sum(last_valid == "")
    if (n_all_unclass > 0) {
      warning(n_all_unclass, " row(s) are fully unclassified up to rank '", taxrank, "'.")
    }
  }

  # Update target column
  taxa[is_unclass, idx] <- ifelse(last_valid == "", "Unclassified",
    paste("Unclassified", last_valid)
  )

  return(taxa)
}


#' Fill all unclassified taxonomic ranks using the last known parent clade
#'
#' This internal function replaces missing or empty taxonomic entries (including `NA`, `""`, `" "`, `"\t"`)
#' with "Unclassified". It then processes all taxonomic ranks sequentially, relabeling "Unclassified"
#' entries as "Unclassified {Last_Known_Parent}" using the last known valid parent clade from higher ranks.
#'
#' The function propagates the most recent non-unclassified taxon name across lower ranks,
#' ensuring consistent hierarchical context for unclassified taxa. Fully unclassified rows
#' (i.e., all ranks are "Unclassified") remain labeled as "Unclassified".
#'
#' @param taxa A `data.frame` or `matrix` where columns represent taxonomic ranks
#'   (e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species) in hierarchical order,
#'   and rows represent taxa.
#' @param silent A `logical` value indicating whether to suppress warnings.
#'   If `FALSE`, a warning is issued if any rows are fully unclassified across all ranks.
#'
#' @return A `matrix` of the same dimensions as input, with:
#'   - All missing/empty values replaced by "Unclassified".
#'   - "Unclassified" entries in each rank relabeled as "Unclassified {Parent}"
#'     where `{Parent}` is the last non-unclassified taxon in the row up to that point.
#'   - Fully unclassified rows remain as "Unclassified" in all columns.
#'   - All original columns preserved (no column subsetting).
#'
#' @examples
#' taxa <- data.frame(
#'   Kingdom = c("Kingdom1", "Kingdom2", "Kingdom3", "Kingdom4", "Kingdom5", "Kingdom6", ""),
#'   Phylum = c("Phylum1", "", "Phylum3", "Phylum4", "Phylum5", "", NA),
#'   Class = c("Class1", "Class2", " ", "\t", NA, "Class6", ""),
#'   Order = c("Order1", "Order2", NA, "Order4", NA, "", " ")
#' )
#' rownames(taxa) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5", "OTU6", "OTU7")
#'
#' .fill_unclassified_all(taxa, silent = FALSE)
#'
#' @keywords internal
#' @noRd
.fill_unclassified_all <- function(taxa, silent = TRUE) {
  stopifnot(is.data.frame(taxa) || is.matrix(taxa))

  # Convert to character matrix for speed
  taxa <- as.matrix(taxa)

  # Replace all NA's or empty strings
  empty_vals <- c(NA, "", " ", "\t")
  taxa[taxa %in% empty_vals] <- "Unclassified"

  # Track last valid (non-unclassified) name per row across ranks
  last_valid <- rep("", nrow(taxa))

  for (j in seq_len(ncol(taxa))) {
    col <- taxa[, j]
    is_unclass <- col == "Unclassified"
    # Replace unclassified with "Unclassified <Parent>" using last valid name
    col[is_unclass] <- ifelse(last_valid[is_unclass] == "",
      "Unclassified",
      paste("Unclassified", last_valid[is_unclass])
    )
    taxa[, j] <- col

    # Update last_valid with names of "true" taxa (not "Unclassified" or
    # "Unclassified <Parent>") so they propagate to lower ranks
    true_name <- !grepl("^Unclassified($| )", taxa[, j])
    last_valid[true_name] <- taxa[true_name, j]
  }

  if (!silent) {
    all_unclass <- apply(taxa, 1, function(x) all(grepl("^Unclassified($| )", x)))
    n_all_unclass <- sum(all_unclass)
    if (n_all_unclass > 0) {
      warning(n_all_unclass, " row(s) are fully unclassified across all ranks.")
    }
  }

  return(taxa)
}


#' Apply Total Sum Scaling (TSS) to numeric columns
#'
#' @param x A `data.frame` or `matrix` containing numeric columns to be scaled.
#' @param presence_threshold `Numeric` value indicating the minimum threshold for
#' presence; values below this are set to 0. Default is 0.
#' @param prop_scale Scaling factor applied after normalization. Default is 100
#' (i.e., convert to percentages).
#'
#' @returns The input object x with numeric columns transformed by TSS.
#' Non-numeric columns remain unchanged.
#' @keywords internal
#' @noRd
.apply_tss <- function(x, presence_threshold = 0, prop_scale = 100) {
  if (any(x < 0, na.rm = TRUE)) {
    warning("Data contains negative entries. Result from applying TSS may not make sense.")
  }

  # Handle data.frame and matrix differently
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("Matrix must be numeric.")
    # # Vectorized thresholding
    # x[x < presence_threshold] <- 0
    # col_sums <- colSums(x, na.rm = TRUE)
    # # Replace zero sums to avoid NaN
    # col_sums[col_sums == 0] <- 1
    # x <- sweep(x, 2, col_sums, "/") * prop_scale
    # x[is.nan(x)] <- 0

    # For speed, convert to data.frame, process, then convert back as data.frames
    # store columns as list elements, allowing fast column access (matrices are
    # vectors with dimensions; column extraction requires more computation)
    x_df <- as.data.frame(x)
    x_df[] <- lapply(x_df, function(col) {
      col_thresh <- ifelse(col >= presence_threshold, col, 0)
      if (all(col_thresh == 0, na.rm = TRUE)) {
        return(rep(0, length(col)))
      }
      col_norm <- col_thresh / sum(col_thresh, na.rm = TRUE) * prop_scale
      ifelse(is.nan(col_norm), 0, col_norm)
    })
    x <- as.matrix(x_df)
  } else {
    # Get all numeric columns
    is_numeric <- vapply(x, is.numeric, logical(1L))

    x[is_numeric] <- lapply(x[is_numeric], function(col) {
      col_thresh <- ifelse(col >= presence_threshold, col, 0)

      # Early exit if all values are zero after thresholding to avoid unnecessary
      # computation
      if (all(col_thresh == 0, na.rm = TRUE)) {
        return(rep(0, length(col_thresh)))
      }

      col_norm <- col_thresh / sum(col_thresh, na.rm = TRUE) * prop_scale
      # Replace NaN values resulting from division (e.g., 0/0) with 0
      ifelse(is.nan(col_norm), 0, col_norm)
    })
  }
  return(x)
}


#' Agglomerate filtered_taxa in microEDA info slot
#'
#' @keywords internal
#' @noRd
.agglomerate_filtered_taxa <- function(data_list,
                                       tax_rank,
                                       rm_missing = FALSE,
                                       transform = c("None", "TSS"),
                                       add_prefix = FALSE) {
  # Input validation
  if (!is.list(data_list) || !all(c("filtered_otu", "filtered_tax") %in% names(data_list))) {
    stop("'data_list' must be a list with 'filtered_otu' and 'filtered_tax' matrices.")
  }

  otu_mat <- data_list$filtered_otu
  tax_mat <- data_list$filtered_tax

  if (!is.matrix(otu_mat) || !is.matrix(tax_mat)) {
    stop("'filtered_otu' and 'filtered_tax' must be matrices.")
  }

  if (nrow(otu_mat) != nrow(tax_mat)) {
    stop("Number of taxa in filtered_otu and filtered_tax must match.")
  }

  if (!identical(rownames(otu_mat), rownames(tax_mat))) {
    stop("Taxa in filtered_otu and filtered_tax must have the same row names and order.")
  }

  if (!tax_rank %in% colnames(tax_mat)) {
    stop("Selected taxonomic rank '", tax_rank, "' not found in filtered_tax.")
  }

  transform <- match.arg(transform, choices = c("None", "TSS"))
  rank_index <- which(colnames(tax_mat) == tax_rank)

  # Handle missing values
  if (rm_missing) {
    empty_vals <- c(NA, "", " ", "\t", "Unassigned", "NA")
    keep_rows <- !(is.na(tax_mat[, rank_index]) | tax_mat[, rank_index] %in% empty_vals)
    tax_mat <- tax_mat[keep_rows, , drop = FALSE]
    otu_mat <- otu_mat[keep_rows, , drop = FALSE]
  } else {
    tax_mat <- .fill_unclassified_all(tax_mat)
  }

  # Add prefixes if requested
  if (add_prefix) {
    prefix_map <- c(
      "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
      "Order" = "o__", "Family" = "f__", "Genus" = "g__",
      "Species" = "s__", "Strain" = "t__"
    )
    for (rank in names(prefix_map)) {
      if (rank %in% colnames(tax_mat)) {
        tax_mat[, rank] <- ifelse(
          grepl("^k__|p__|c__|o__|f__|g__|s__|t__", tax_mat[, rank]),
          tax_mat[, rank],
          paste0(prefix_map[rank], tax_mat[, rank])
        )
      }
    }
  }

  # Create group ID based on full lineage up to tax_rank
  tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
  tax_df$tax_ID <- apply(
    tax_df[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  # Agglomerate abundances
  aggregated <- rowsum(as.matrix(otu_mat), group = tax_df$tax_ID, reorder = FALSE)

  # Apply TSS if requested
  if (transform == "TSS") {
    aggregated <- .apply_tss(aggregated, prop_scale = 100)
  }

  # Update tax_table: keep one row per tax_ID

  # Extract full taxonomy up to the target rank
  tax_subset <- tax_df[, 1:rank_index, drop = FALSE]

  # To match with aggregated tax_IDs (the grouped names at target rank) and to
  # ensure row order matches aggregated
  tax_subset$tax_ID <- apply(
    tax_subset[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  new_tax_df <- tax_subset[!duplicated(tax_subset$tax_ID), , drop = FALSE]
  new_tax_df <- new_tax_df[match(rownames(aggregated), new_tax_df$tax_ID), ]
  new_tax_mat <- as.matrix(new_tax_df[, -ncol(new_tax_df)]) # Since tax_ID is last column
  rownames(aggregated) <- rownames(new_tax_mat)

  # Return updated list
  return(list(filtered_otu = aggregated, filtered_tax = new_tax_mat))
}


#' Agglomerate Taxa by Taxonomic Rank
#'
#' Collapses a microbiome dataset (`microEDA` or `phyloseq` object) by aggregating
#' features of the same taxonomy into higher-level taxonomic groups based on a
#' specified rank.
#' Optionally removes unclassified taxa, applies total sum scaling (TSS), or
#' adds standard taxonomic rank prefixes (e.g., `p__`).
#'
#' @param me A `microEDA` or `phyloseq` object containing abundance and taxonomic data.
#' @param tax_rank `Character` string specifying the taxonomic rank for agglomeration
#'   (e.g., "Genus", "Family"). Must exist in the tax_table of `me`.
#' @param rm_missing `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values using lineage context.
#' @param transform `Character.` Transformation to apply to abundances after agglomeration.
#'   One of `"None"` (no transformation) or `"TSS"` (Total Sum Scaling to relative abundance).
#' @param add_prefix `Logical`. If `TRUE`, adds QIIME-style prefixes (e.g., `k__`, `p__`)
#'   to taxonomic labels.
#'
#' @return Returns an object of the same class as input (`microEDA` or `phyloseq`),
#'   with taxa agglomerated at the specified rank. Preserves sample data, phylogenetic
#'   tree (pruned), and reference sequences if present in the input.
#'
#' @details
#' This function performs safe agglomeration by grouping taxa based on their full
#' lineage up to the target rank, reducing the risk of incorrectly merging distinct clades.
#' A warning is issued if multiple distinct higher-rank lineages map to the same group
#' at the target level, indicating potential annotation inconsistencies.
#'
#' @note This implementation is significantly faster than [`phyloseq::tax_glom`],
#' especially on large datasets, due to efficient use of vectorized operations.
#' Benchmarks show speedups of over 100x compared to `tax_glom`.
#'
#' @examples
#' \dontrun{
#' # Example with a phyloseq object
#' data("GlobalPatterns", package = "phyloseq")
#' agglom <- agglomerate_taxa(GlobalPatterns,
#'   tax_rank = "Family",
#'   rm_missing = TRUE, transform = "TSS",
#'   add_prefix = TRUE
#' )
#' print(agglom)
#' }
#' @export
agglomerate_taxa <- function(me,
                             tax_rank,
                             rm_missing = FALSE,
                             transform = c("None", "TSS"),
                             add_prefix = FALSE) {
  # Input validation
  if (!inherits(me, "microEDA") && !inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }

  transform <- match.arg(transform, choices = c("None", "TSS"))

  # Extract data using phyloseq accessors
  abund_table <- phyloseq::otu_table(me)
  tax_tab <- phyloseq::tax_table(me)

  if (is.null(tax_tab)) {
    stop('The "tax_table" slot must not be empty.')
  }

  # Check consistency
  if (nrow(abund_table) != nrow(tax_tab)) stop("Abundance and taxonomic table do not have the same number of taxa!")

  if (!tax_rank %in% colnames(tax_tab)) stop("Selected taxonomic rank '", tax_rank, "' not found in tax_table.")

  # Convert to data frames for easier manipulation
  tax_df <- as.data.frame(as(tax_tab, "matrix"), stringsAsFactors = FALSE)
  abund_df <- as.data.frame(as(abund_table, "matrix"), stringsAsFactors = FALSE)
  if (!identical(rownames(abund_df), rownames(tax_df))) stop("Taxa in abundance and taxonomic table are not in the same order!")

  rank_index <- which(colnames(tax_df) == tax_rank)

  # Handle unclassified/missing taxa
  if (rm_missing) {
    empty_vals <- c(NA, "", " ", "\t", "Unassigned", "NA")
    keep_rows <- !(is.na(tax_df[, rank_index]) | tax_df[, rank_index] %in% empty_vals)
    tax_df <- tax_df[keep_rows, , drop = FALSE]
    abund_df <- abund_df[keep_rows, , drop = FALSE]
  } else {
    # Fill missing taxa with "Unclassified" <Parent>
    tax_df <- as.data.frame(.fill_unclassified_all(tax_df), stringsAsFactors = FALSE)
  }

  # Apply prefixes to all taxonomic ranks if add_prefix = TRUE
  if (add_prefix) {
    # Define prefix mapping
    prefix_map <- c(
      "Kingdom" = "k__", "Phylum" = "p__", "Class" = "c__",
      "Order" = "o__", "Family" = "f__", "Genus" = "g__",
      "Species" = "s__", "Strain" = "t__"
    )
    for (rank in names(prefix_map)) {
      if (rank %in% colnames(tax_df)) {
        # Avoid double-prefixing
        tax_df[[rank]] <- ifelse(
          grepl("^k__|p__|c__|o__|f__|g__|s__|t__", tax_df[[rank]]),
          tax_df[[rank]], # Keep if already prefixed
          paste0(prefix_map[rank], tax_df[[rank]])
        )
      }
    }
  }

  # Add tax_ID for grouping
  tax_df$tax_ID <- tax_df[[tax_rank]]

  # Check for conflicting higher-level taxonomy
  tax_check <- tax_df |>
    dplyr::group_by(tax_ID) |>
    dplyr::summarise(
      across(1:all_of(tax_rank),
        ~ dplyr::n_distinct(.x),
        .names = "n_{col}"
      ),
      .groups = "drop"
    )

  # Check if any rank has more than one distinct value per tax_ID
  conflicts <- dplyr::filter(tax_check, if_any(starts_with("n_"), ~ .x > 1))

  if (nrow(conflicts) > 0) {
    warning(
      "Conflicting taxonomy detected for ",
      nrow(conflicts),
      " tax_ID(s). Some taxa grouped at '",
      tax_rank,
      "' have inconsistent higher-level classifications. These are:\n",
      toString(conflicts$tax_ID)
    )
  }

  # Collapse each lineage into a single string to serve as groups
  tax_df$tax_ID <- apply(
    tax_df[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )
  # Ensure tax_ID is same for tax and abund tables
  # abund_df$tax_ID <- tax_df$tax_ID

  group_ids <- tax_df$tax_ID
  # names(group_ids) <- rownames(tax_df)

  # Agglomerate abundances by tax_ID
  aggregated <- rowsum(as(abund_df, "matrix"), group = group_ids, reorder = FALSE)
  # aggregated <- abund_df |>
  #   # tibble::rownames_to_column("OTU") |>
  #   dplyr::group_by(tax_ID) |>
  #   dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  # Apply TSS (Total Sum Scaling / Relative Abundance) if requested
  if (transform == "TSS") {
    aggregated <- .apply_tss(aggregated, prop_scale = 100)
  }

  # Extract full taxonomy up to the target rank
  tax_subset <- tax_df[, 1:rank_index, drop = FALSE]

  # To match with aggregated tax_IDs (the grouped names at target rank) and to
  # ensure row order matches aggregated$tax_ID
  # tax_subset$tax_ID <- tax_subset[[tax_rank]]
  tax_subset$tax_ID <- apply(
    tax_subset[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  # Aggregate: keep only one row per tax_ID, keeping the first occurrence (since
  # all rows with same tax_ID at target rank should have same higher ranks, ideally)
  new_tax_df <- tax_subset[!duplicated(tax_subset$tax_ID), , drop = FALSE]

  # Reorder to match the order in 'aggregated'
  # new_tax_df <- new_tax_df[match(aggregated$tax_ID, new_tax_df$tax_ID), ]
  new_tax_df <- new_tax_df[match(rownames(aggregated), new_tax_df$tax_ID), ]

  # Remove tax_ID column and convert to matrix for tax_table
  new_tax_df <- new_tax_df[, -which(names(new_tax_df) == "tax_ID"), drop = FALSE]
  tax_tbl <- phyloseq::tax_table(as.matrix(new_tax_df))

  # Update row names and remove tax_ID for OTU table
  rownames(aggregated) <- rownames(tax_tbl)
  # aggregated <- aggregated[, -which(names(aggregated) == "tax_ID")]
  otu_tbl <- phyloseq::otu_table(as.matrix(aggregated), taxa_are_rows = TRUE)

  # Construct new phyloseq object
  phy_new <- phyloseq::phyloseq(otu_tbl, tax_tbl)

  # Preserve sample data if present
  if (!is.null(phyloseq::sample_data(me, errorIfNULL = FALSE))) {
    phy_new <- phyloseq::merge_phyloseq(phy_new, sample_data(me))
  }

  # Preserve and prune phylogenetic tree if present
  if (!is.null(phyloseq::phy_tree(me, errorIfNULL = FALSE))) {
    tree <- phyloseq::phy_tree(me)
    tree <- ape::drop.tip(tree, tip = setdiff(tree$tip.label, phyloseq::taxa_names(phy_new)))
    phy_new <- phyloseq::merge_phyloseq(phy_new, tree)
  }

  # Preserve refseq if present
  if (!is.null(phyloseq::refseq(me, errorIfNULL = FALSE))) {
    refseq_data <- phyloseq::refseq(me)
    refseq_data <- refseq_data[phyloseq::taxa_names(phy_new)] # prune
    phy_new <- phyloseq::merge_phyloseq(phy_new, refseq_data)
  }

  if (inherits(me, "microEDA")) {
    # Reconstruct microEDA
    agglomerated_obj <- new("microEDA", phy_new)

    # Preserve info slot
    info(agglomerated_obj) <- info(me)
    # Update taxrank
    taxrank(agglomerated_obj) <- tax_rank

    # Update transform if necessary
    if (!transform == "None") transforms(agglomerated_obj) <- transform

    # If present, also agglomerate filtered taxa
    other_tax <- filtered_taxa(me)

    if (!is.null(other_tax)) {
      agg_filtered <- .agglomerate_filtered_taxa(other_tax,
        tax_rank = tax_rank,
        rm_missing = rm_missing,
        transform = transform,
        add_prefix = add_prefix
      )
      filtered_taxa(agglomerated_obj) <- agg_filtered
    }
  } else {
    agglomerated_obj <- phy_new
  }
  return(agglomerated_obj)
}
