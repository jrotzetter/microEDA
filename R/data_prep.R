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


#' Apply Data Transformation to Abundance Matrix
#'
#' Internal function to apply normalisations or transformations such as
#' Total Sum Scaling (TSS) to the OTU table of a phyloseq-like object. Optionally
#' integrates and transforms an auxiliary matrix of filtered taxa if present.
#'
#' The function combines the main OTU table with the filtered taxa matrix
#' (if available in a `microEDA` object) by aligning on shared samples. It then
#' applies the specified transformation (e.g., TSS for relative abundance). After
#' transformation, the matrix is split back into original components, row names
#' are restored, and the transformed matrices are reassigned. The transformation
#' type is stored in a microEDA's metadata.
#' @param me A `microEDA` or `phyloseq` object.
#' @param transform `Character` string specifying the transformation. One of:
#'  `"None"`, `"TSS"`.
#'
#' @return The modified object `me` with transformed abundance matrices and updated
#'  transformation metadata.
#' @details
#' - TSS Transformation: Scales counts per sample to sum to a constant
#'  (e.g., relative abundance x 100). Implemented via internal `.apply_tss` function.
#' - Filtered Taxa Integration: Only performed if `inherits(me, "microEDA")` and
#' `filtered_taxa(me)$filtered_otu` exists. Matrices are merged row-wise,
#'  transformed, then split back using stored dimensions.
#' - Column Alignment: Only samples present in both matrices are used. An error
#'  is thrown if no common columns exist.
#'
#' @keywords internal
#' @noRd
.apply_transformation <- function(me, transform = c("None", "TSS")) {
  transform <- match.arg(transform, choices = c("None", "TSS"))

  abund_mat <- as(phyloseq::otu_table(me), "matrix")

  abund_rows <- rownames(abund_mat)
  n_abund <- nrow(abund_mat)

  if (inherits(me, "microEDA")) {
    other_mat <- filtered_taxa(me)$filtered_otu
    if (!is.null(other_mat)) {
      # Find common columns
      common_cols <- intersect(colnames(abund_mat), colnames(other_mat))

      if (length(common_cols) == 0) {
        stop("No common columns between 'otu_table' and 'filtered_taxa'.")
      }

      # Subset both matrices to shared columns
      abund_mat <- abund_mat[, common_cols, drop = FALSE]
      other_mat <- other_mat[, common_cols, drop = FALSE]

      other_rows <- rownames(other_mat)
      n_other <- nrow(other_mat)

      # Combine
      abund_mat <- rbind(abund_mat, other_mat)
    }
  }
  # Apply transformation
  switch(transform,
    None = return(me),
    TSS = {
      abund_mat <- .apply_tss(abund_mat, presence_threshold = 0, prop_scale = 100)
    }
  )

  # Split back
  if (inherits(me, "microEDA") && !is.null(other_mat)) {
    abund_mat_out <- abund_mat[seq_len(n_abund), , drop = FALSE]
    other_mat_out <- abund_mat[seq_len(n_other) + n_abund, , drop = FALSE]

    # Restore row names
    rownames(abund_mat_out) <- abund_rows
    rownames(other_mat_out) <- other_rows

    filtered_taxa(me)$filtered_otu <- other_mat_out
  } else {
    abund_mat_out <- abund_mat
  }

  otu_table(me) <- otu_table(abund_mat_out, taxa_are_rows = TRUE)
  if (inherits(me, "microEDA")) transforms(me) <- transform

  return(me)
}


#' Agglomerate filtered_taxa in microEDA info slot
#'
#' @keywords internal
#' @noRd
.agglomerate_filtered_taxa <- function(data_list,
                                       tax_rank,
                                       rm_missing = FALSE,
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
    tax_mat <- .add_taxonomy_prefix(tax_mat)
  }

  # Create group ID based on full lineage up to tax_rank
  tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
  tax_df$tax_ID <- apply(
    tax_df[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  # Agglomerate abundances by tax_ID
  aggregated <- rowsum(as.matrix(otu_mat), group = tax_df$tax_ID, reorder = FALSE)

  # Update tax_table: keep one row per tax_ID

  # Extract full taxonomy up to the target rank
  tax_subset <- tax_df[, 1:rank_index, drop = FALSE]

  # To match with aggregated tax_IDs (the grouped names at target rank) and to
  # ensure row order matches aggregated
  tax_subset$tax_ID <- apply(
    tax_subset[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  # Add total abundance over all samples for each OTU
  abund_sums <- rowSums(otu_mat)

  # Get rownames (ASVs/OTUs) with highest abundance within a tax_ID, which
  # should be biologically more representative than using the first occurring
  # OTU of a tax_ID
  max_rownames <- tapply(seq_along(tax_subset$tax_ID), tax_subset$tax_ID, function(i) i[which.max(abund_sums[i])])
  max_rownames <- unname(unlist(max_rownames))

  # Then use these to keep only one row per tax_ID
  new_tax_df <- tax_subset[max_rownames, ]

  # Reorder to match the order in 'aggregated'
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
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae").
#' @param transform `Character.` Transformation to apply to abundances after agglomeration.
#'   One of `"None"` (no transformation) or `"TSS"` (Total Sum Scaling to relative abundance).
#'   If `filtered_taxa` are present in a `microEDA` object, they will be included in the transformation calculation.
#' @param add_prefix `Logical`. If `TRUE`, adds QIIME-style prefixes (e.g., `k__`, `p__`)
#'   to taxonomic labels.
#'
#' @return Returns an object of the same class as input (`microEDA` or `phyloseq`),
#'   with taxa agglomerated at the specified rank. Preserves sample data, phylogenetic
#'   tree (pruned), and reference sequences if present in `me`.
#'
#' @details
#' This function performs safe agglomeration by grouping taxa based on their full
#' lineage up to the target rank, reducing the risk of incorrectly merging distinct clades.
#' A warning is issued if multiple distinct higher-rank lineages map to the same group
#' at the target level, indicating potential annotation inconsistencies.
#'
#' When multiple sequences (ASVs/OTUs) belong to the same taxonomic group at the
#' target rank, the ID assigned to the agglomerated feature is taken from
#' the most abundant sequence (summed across samples) within that group. This
#' should ensure that the dominant biological variant drives ASV/OTU labeling,
#' enhancing representativeness.
#'
#' @note This implementation is significantly faster than [`phyloseq::tax_glom`],
#' especially on large datasets, due to efficient use of vectorized operations.
#' Benchmarks show speedups of 10x to over 100x compared to `tax_glom`.
#'
#' @examples
#' # Example with a phyloseq object
#' data("GlobalPatterns", package = "phyloseq")
#' agglom <- agglomerate_taxa(GlobalPatterns,
#'   tax_rank = "Phylum",
#'   rm_missing = TRUE, transform = "TSS",
#'   add_prefix = TRUE
#' )
#' agglom
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

  tax_rank <- .get_full_tax_rank(tax_rank) # In case it was abbreviated

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
    tax_df <- .add_taxonomy_prefix(tax_df)
  }

  # Check for conflicting higher-level taxonomy
  suppressMessages(check_taxonomic_consistency(tax_df, tax_rank = tax_rank, detailed_report = FALSE))

  # Add tax_ID for grouping
  # Collapse each lineage into a single string to serve as groups
  tax_df$tax_ID <- apply(
    tax_df[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  group_ids <- tax_df$tax_ID
  # Agglomerate abundances by tax_ID
  aggregated <- rowsum(as(abund_df, "matrix"), group = group_ids, reorder = FALSE)

  # Extract full taxonomy up to the target rank
  tax_subset <- tax_df[, 1:rank_index, drop = FALSE]

  # To match with aggregated tax_IDs (the grouped names at target rank) and to
  # ensure row order matches aggregated$tax_ID
  tax_subset$tax_ID <- apply(
    tax_subset[, 1:rank_index, drop = FALSE], 1,
    function(x) paste(x, collapse = ";")
  )

  # Add total abundance over all samples for each OTU
  abund_sums <- rowSums(abund_df)

  # Get rownames (ASVs/OTUs) with highest abundance within a tax_ID, which
  # should be biologically more representative than using the first occurring
  # OTU of a tax_ID
  max_rownames <- tapply(seq_along(tax_subset$tax_ID), tax_subset$tax_ID, function(i) i[which.max(abund_sums[i])])
  max_rownames <- unname(unlist(max_rownames))

  # Then use these to keep only one row per tax_ID
  new_tax_df <- tax_subset[max_rownames, ]

  # Reorder to match the order in 'aggregated'
  new_tax_df <- new_tax_df[match(rownames(aggregated), new_tax_df$tax_ID), ]

  # Remove tax_ID column and convert to matrix for tax_table
  new_tax_df <- new_tax_df[, -which(names(new_tax_df) == "tax_ID"), drop = FALSE]
  tax_tbl <- phyloseq::tax_table(as.matrix(new_tax_df))

  # Update row names and remove tax_ID for OTU table
  rownames(aggregated) <- rownames(tax_tbl)
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

    # If present, also agglomerate filtered taxa
    other_tax <- filtered_taxa(me)

    if (!is.null(other_tax)) {
      agg_filtered <- .agglomerate_filtered_taxa(other_tax,
        tax_rank = tax_rank,
        rm_missing = rm_missing,
        add_prefix = add_prefix
      )
      filtered_taxa(agglomerated_obj) <- agg_filtered
    }
  } else {
    agglomerated_obj <- phy_new
  }

  # Apply normalisation / transformation if requested
  if (!transform == "None") agglomerated_obj <- .apply_transformation(agglomerated_obj, transform = transform)

  return(agglomerated_obj)
}


#' Prepare taxonomic profile for plotting
#'
#' Aggregates and filters taxa to a manageable number for visualization,
#' applying abundance and prevalence filters, rank-based agglomeration,
#' and iterative thresholding to limit the number of displayed taxa.
#'
#' @returns A long-format data frame with columns:
#'  \itemize{
#'    \item \code{Taxon}: Taxon name, with "Other" labeled by threshold.
#'    \item \code{Sample}: Sample identifier.
#'    \item \code{Abundance}: Abundance value.
#'  }
#'
#' @details
#' Note that the abundance of the "Other" category in the final output may
#' exceed the effective threshold because it represents the sum of multiple
#' individual taxa, each of which had relative abundances below the threshold
#' after initial filtering and agglomeration.
#' @keywords internal
#' @noRd
.prepare_tax_profile <- function(me,
                                 tax_rank,
                                 min_abundance = 0,
                                 min_prevalence = 0,
                                 ntaxa = 30,
                                 group_var = NULL,
                                 abundance_criterion = c("prevalence", "mean"),
                                 filter_by_group = FALSE,
                                 group_requirement = c("any", "all"),
                                 keep_filtered = TRUE,
                                 rm_missing = FALSE,
                                 transform = c("None", "TSS"),
                                 add_prefix = FALSE,
                                 process_taxon = TRUE) {
  # Input validation
  if (!inherits(me, "microEDA") && !inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }

  abundance_criterion <- match.arg(abundance_criterion, choices = c("prevalence", "mean"))
  group_requirement <- match.arg(group_requirement, choices = c("any", "all"))
  transform <- match.arg(transform, choices = c("None", "TSS"))

  if (!tax_rank %in% phyloseq::rank_names(me)) {
    stop(
      "Rank '", tax_rank, "' not present in tax_table. Available ranks: ",
      paste(phyloseq::rank_names(tax_rank), collapse = ", ")
    )
  }

  # Check for nonsensical filtering conditions
  if ((min_abundance == 0) != (min_prevalence == 0)) {
    stop("One of 'min_abundance' or 'min_prevalence' is zero while the other is not. This may lead to unintended filtering behavior as the zero threshold will pass all values.")
  }

  if (any(phyloseq::otu_table(me) < 0, na.rm = TRUE)) {
    stop("Abundance data contains negative values.")
  }

  group_var_arg <- if (filter_by_group) group_var else NULL

  me <- filter_features(me,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    group_var = group_var_arg,
    abundance_criterion = abundance_criterion,
    group_requirement = group_requirement,
    keep_filtered = keep_filtered
  )

  me <- agglomerate_taxa(me,
    tax_rank = tax_rank,
    rm_missing = rm_missing,
    transform = transform,
    add_prefix = add_prefix
  )

  # Extract OTU table, then append tax_ID column at target rank
  abund_tab <- as.data.frame(phyloseq::otu_table(me))
  tax_ids <- as.data.frame(phyloseq::tax_table(me), stringsAsFactors = FALSE)[[tax_rank]]
  abund_tab <- data.frame(tax_ID = tax_ids, abund_tab, stringsAsFactors = FALSE)

  # Get filtered out taxa for "Other" row (only available for microEDA)
  if (inherits(me, "microEDA") && !is.null(filtered_taxa(me))) {
    other_row <- .collapse_other(me)$other_row
    other_row <- data.frame(tax_ID = "Other", other_row, stringsAsFactors = FALSE)

    # Now append "Other" row to OTU table
    abund_tab <- rbind(abund_tab, other_row)
  }

  # Check if data is relative abundances or counts for plot labeling
  if (suppressWarnings(.is_proportion(abund_tab))) {
    is_rel_abund <- TRUE
  } else if (suppressWarnings(.is_counts(abund_tab))) {
    is_rel_abund <- FALSE
  } else {
    stop("Data is neither counts nor relative abundance - data was likely transformed.")
  }

  # Pass to incrementing filter to reduce ntaxa
  abund_tab_reduced <- .apply_incrementing_filter(abund_tab,
    ntaxa = ntaxa,
    initial_threshold = 0
  )

  tax_abund <- abund_tab_reduced$data |>
    dplyr::rename(Taxon = "tax_ID") |>
    # Sum abundances of duplicate taxa across higher-level lineage variants.
    # This is necessary for taxa with conflicting higher-level inconsistencies
    # due to several entries for the same taxon existing per sample, introducing
    # bias to summary statistics for the taxon (mean abundance, standard deviation, prevalence)
    dplyr::group_by(.data$Taxon) |>
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop") |>
    tidyr::pivot_longer(-.data$Taxon, names_to = "Sample", values_to = ifelse(is_rel_abund, "Relative_Abundance", "Abundance")) |>
    dplyr::arrange(.data$Sample)

  # Process tax_ID for plotting (markdown)
  if (process_taxon) {
    tax_abund <- tax_abund |>
      dplyr::mutate(
        Taxon = stringi::stri_replace_first_regex(
          .data$Taxon,
          "(.*)_unclassified", "Unclassified *$1*"
        ),
        Taxon = stringi::stri_replace_first_regex(
          .data$Taxon,
          "Unclassified (.*)", "Unclassified *$1*"
        ),
        Taxon = stringi::stri_replace_first_regex(
          .data$Taxon,
          "^(\\S*)$", "*$1*"
        ),
        Taxon = stringi::stri_replace_all_regex(.data$Taxon, "_", " "),
        Taxon = dplyr::if_else(.data$Taxon == "*Other*", "Other", .data$Taxon),
        Taxon = dplyr::if_else(.data$Taxon == "*UNCLASSIFIED*", "UNCLASSIFIED", .data$Taxon)
      )
  }
  return(tax_abund)
}


#' Recode values in a data frame column using named replacements
#'
#' This internal function renames values in a specified column of a data frame
#' based on a named list of replacements. It uses `dplyr::recode()` with splicing
#' (`!!!`) to apply the mapping.
#'
#' @param df A `data.frame`.
#' @param column A `character` string specifying the name of the column to modify.
#' @param values A named `list` or named `vector` where names are old values and values are new values.
#'
#' @return The input data frame with the specified column's values renamed according to `values`.
#'
#' @examples
#' df <- data.frame(x = c("a", "b", "c"))
#' new_vals <- list(a = "Alfa", b = "Bravo", c = "Charlie")
#' .rename_values(df, "x", new_vals)
#'
#' @keywords internal
#' @noRd
.rename_values <- function(df, column, values) {
  col <- df[[column]]
  col_new <- dplyr::recode(col, !!!values)
  df[[column]] <- col_new
  return(df)
}


#' @keywords internal
#' @noRd
.specify_decimal <- function(x, k = 3) {
  # will return a character string, only for printing/plotting purposes
  trimws(format(round(x, k), nsmall = k))
}
