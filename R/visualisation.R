#' Generate a Named Color Vector for Categorical Levels
#'
#' This internal function creates a color vector with one color per unique level
#' of a specified categorical variable in a data frame. Currently, colors are
#' drawn from a RColorBrewer palette via `colorRampPalette`. Optionally, the
#' color for a category named "Other (<...)" can be set to "darkgrey".
#'
#' @param df A `data.frame` containing the categorical variable.
#' @param level Unquoted column name (tidy evaluation) representing the
#'              categorical variable for which colors are generated.
#' @param other_cat `Logical`. If `TRUE`, sets the color of any category whose
#'                  name matches the pattern "Other (<...)" to "darkgrey".
#' @param pal `Character` string specifying the RColorBrewer palette to use
#'            (default is "Paired"). Must be a valid palette name with at least
#'            12 colors.
#'
#' @return A named character vector of color codes, where:
#'         - Names are the unique values of `df[[level]]`.
#'         - Values are hex color codes interpolated from the chosen palette.
#'         - Entries with `NA` names are removed.
#'         - If `other_cat = TRUE`, the "Other" category is colored "darkgrey".
#'
#' @details
#' - The function uses `dplyr::ungroup()` to ensure grouping does not affect
#'   distinct level extraction.
#' - `{{ level }}` uses tidy evaluation (from `rlang`) to allow unquoted column names.
#' - `grDevices::colorRampPalette` interpolates colors if more levels exist than
#'   colors in the base palette (e.g., more than 12).
#' - Categories with `NA` as a name are filtered out before return.
#'
#' @importFrom dplyr ungroup distinct pull
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @keywords internal
#' @noRd
.get_col_vector <- function(df, level, other_cat = FALSE, pal = "Paired") {
  num_colors <- df |>
    dplyr::ungroup() |>
    dplyr::distinct({{ level }}) |>
    nrow()

  col_vector <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, pal))(num_colors)

  names(col_vector) <- df |>
    dplyr::ungroup() |>
    dplyr::distinct({{ level }}) |>
    dplyr::pull()

  if (other_cat) {
    index <- grep("^Other\\s\\(<.*\\)$", names(col_vector))
    col_vector[index] <- "darkgrey"
  }
  col_vector <- col_vector[!is.na(names(col_vector))]
  return(col_vector)
}


#' Plot Taxonomic Composition of a microEDA/phyloseq Object
#'
#' Creates a ggplot2-based stacked bar chart showing the abundance of taxa
#' across samples, optionally grouped and faceted by metadata variables.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU table, taxonomy, and sample data.
#' @param tax_rank `Character` string specifying the taxonomic rank to plot
#'                 (e.g., "Phylum", "Family"). Must be a valid rank in the
#'                 taxonomy table.
#' @param group_var `Character` string (optional) indicating a sample metadata
#'                  variable to group samples.
#' @param facet_var `Character` string (optional) indicating a metadata variable
#'                  for additional faceting/grouping of `group_var`.
#' @param plot_title `Character` string for the main plot title.
#' @param group_labels Named character `vector` mapping old group names to
#'                 new labels (e.g., `c("Old" = "New")`).
#' @param show_samples `Logical`. Control sample label visibility.
#' @param min_abundance `Numeric` value. Minimum abundance threshold for a feature to be retained.
#'   Must be non-negative. Features with abundance below this are considered absent.
#' @param min_prevalence `Numeric` value. Minimum prevalence required for retention.
#'   If value is < 1, interpreted as proportion of samples; otherwise, as absolute number of samples.
#' @param as_relative `Logical`. If `TRUE`, applies relative abundance transformation (TSS)
#'   to the input counts. If `FALSE`, uses raw counts. (Default: `TRUE`).
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in plotting without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param ... Additional arguments for fine-tuning. Can include:
#'   \itemize{
#'     \item `sample_text`: Font size for sample labels (Default: 8).
#'     \item `strip_angle`: Angle of facet strip text (Default: 90).
#'     \item `legend_text`: Font size for legend (Default: 9).
#'     \item `legend_key`: Size of legend keys (Default: 8).
#'     \item `bar_width`: Width of bars (Default: 0.95).
#'     \item `ntaxa`: Number of top taxa to display (Default: 30).
#'     \item `pal`: Color palette name (passed to `RColorBrewer`, default: "Paired").
#'             Must be a valid palette name with at least 12 colors.
#'     \item `abundance_criterion`: `Character` string. Criterion to use for filtering:
#'   \describe{
#'     \item{\code{prevalence}:}{Retain features present in at least `min_prevalence` samples
#'       (within group if `group_var` is used) and with abundance >= `min_abundance`
#'       in those samples.}
#'     \item{\code{mean}:}{Also requires that the mean abundance across samples (or group)
#'       is >= `min_abundance`.}
#'   }
#'   Default: `"prevalence"`.
#'     \item `group_requirement`: `Character` string. When `group_var` is specified, determines
#'   whether the filter criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#'     \item `keep_filtered`: Whether to keep filtered out taxa as "Other"  (Default: `TRUE`). Takes effect only for `microEDA` objects.
#'     \item `rm_missing`: `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae"). (Default: `FALSE`).
#'     \item `process_taxon`: Whether to clean taxon names (remove underscores, return in italic). (Default: `TRUE`).
#'   }
#'
#' @return A `ggplot` object.
#'
#' @details
#' - The function prepares the taxonomic profile by:
#'   \itemize{
#'     \item Filtering taxa by `min_abundance`, `min_prevalence`.
#'     \item Aggregating counts by `tax_rank`.
#'     \item Reducing the total number of taxa to no more than `ntaxa`.
#'     \item Optionally applies relative abundance transformation (TSS) to the data.
#'   }
#' - If `group_var` is specified, samples are grouped by that variable.
#' - If `facet_var` is provided, an additional `facet_wrap` layer is added with markdown formatting.
#' - Sample counts per group are appended to group labels (e.g., "**Control (n=10)**").
#'
#' @examples
#' # Basic plot at Species level
#' mpa <- microEDA(merged_metaphlan_profiles)
#' plot_taxa_barchart(mpa, "Species")
#'
#' @importFrom phyloseq phyloseq sample_data
#' @importFrom dplyr inner_join mutate group_by n_distinct
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace str_replace_all
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual scale_y_continuous facet_grid facet_wrap labs guides
#' @importFrom ggtext element_markdown
#' @importFrom rlang sym
#' @importFrom utils modifyList
#' @importFrom glue glue
#'
#' @export
plot_taxa_barchart <- function(me,
                               tax_rank,
                               group_var = NULL,
                               facet_var = NULL,
                               plot_title = NULL,
                               group_labels = NULL,
                               show_samples = TRUE,
                               min_abundance = 0,
                               min_prevalence = 0,
                               as_relative = TRUE,
                               filter_by_group = FALSE,
                               ...) {
  # Input validation
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

  if (as_relative) {
    transform <- "TSS"
  } else {
    transform <- "None"
  }

  defaults <- list(
    sample_text = 8,
    strip_angle = 90,
    legend_text = 9,
    legend_key = 8,
    bar_width = 0.95,
    ntaxa = 30,
    pal = "Paired",
    abundance_criterion = "prevalence",
    group_requirement = "any", # Only takes effect when group_var != NULL
    keep_filtered = TRUE,
    rm_missing = FALSE,
    process_taxon = TRUE
  )
  # Handle optional ellipsis arguments
  arglist <- list(...)

  arglist <- .warn_invalid_args(allowed = names(defaults), arglist = arglist)

  arglist <- utils::modifyList(defaults, arglist)

  metadata <- NULL
  if (!is.null(group_var)) {
    if (.check_sample_data(me)) {
      if (!(group_var %in% names(phyloseq::sample_data(me)))) {
        stop("'group_var' not found in sample metadata.")
      }

      if (!is.null(group_labels)) {
        phyloseq::sample_data(me) <- .rename_values(phyloseq::sample_data(me), group_var, group_labels)
      }

      metadata <- data.frame(phyloseq::sample_data(me),
        stringsAsFactors = FALSE
      )[, c(group_var, facet_var), drop = FALSE] |>
        .check_var_names() |>
        tibble::rownames_to_column("Sample")

      if (is.factor(metadata[[group_var]])) {
        metadata[[group_var]] <- as.character(metadata[[group_var]])
      }
      group_var_sym <- rlang::sym(group_var)
    } else {
      stop("Can't use 'group_var' without sample metadata!")
    }
  }

  tax_abund <- .prepare_tax_profile(me,
    tax_rank = tax_rank,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    ntaxa = arglist$ntaxa,
    group_var = group_var,
    abundance_criterion = arglist$abundance_criterion,
    filter_by_group = filter_by_group,
    group_requirement = arglist$group_requirement,
    keep_filtered = arglist$keep_filtered,
    rm_missing = arglist$rm_missing,
    transform = transform,
    add_prefix = FALSE,
    process_taxon = arglist$process_taxon,
    calculate_prevalence = FALSE
  )
  col_vector <- .get_col_vector(tax_abund, .data$Taxon, other_cat = TRUE, pal = arglist$pal)

  is_rel_abund <- "Relative_Abundance" %in% colnames(tax_abund)

  y_lab <- ifelse(is_rel_abund, "Relative Abundance (%)", "Abundance")

  if (is_rel_abund) {
    tax_abund <- dplyr::rename(tax_abund, "Abundance" = "Relative_Abundance")
  }

  message("Plotting data...")
  if (is.null(group_var)) {
    ggplot2::ggplot(tax_abund, ggplot2::aes(x = .data$Sample, y = .data$Abundance, fill = .data$Taxon)) +
      ggplot2::geom_col(width = arglist$bar_width, color = "black") +
      ggplot2::scale_fill_manual(tax_rank, values = col_vector) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = NULL,
        y = y_lab,
        title = plot_title
      ) +
      theme_microEDA(
        graph_type = "bar",
        sample_text = arglist$sample_text,
        strip_angle = arglist$strip_angle,
        legend_text = arglist$legend_text,
        legend_key = arglist$legend_key,
        show_samples = show_samples
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
  } else {
    tax_abund <- tax_abund |>
      dplyr::inner_join(metadata, by = "Sample")

    if (!is.null(facet_var)) {
      facet_var_sym <- rlang::sym(facet_var)

      tax_abund <- tax_abund |>
        dplyr::mutate(
          !!group_var_sym := tidyr::replace_na(.data[[group_var]], "Unknown"),
          !!facet_var_sym := stringr::str_replace(.data[[facet_var_sym]], "^(.*)$", "**\\1**")
        ) |>
        dplyr::group_by(.data[[facet_var]], .data[[group_var]]) |>
        dplyr::mutate(num_samples = dplyr::n_distinct(.data$Sample))
    } else {
      tax_abund <- tax_abund |>
        dplyr::mutate(!!group_var_sym := tidyr::replace_na(.data[[group_var]], "Unknown")) |>
        dplyr::group_by(.data[[group_var]]) |>
        dplyr::mutate(num_samples = dplyr::n_distinct(.data$Sample))
    }

    tax_abund <- tax_abund |>
      dplyr::mutate(
        !!group_var_sym := glue::glue("**{.data[[group_var]]}** (n={.data$num_samples})"),
        !!group_var_sym := stringr::str_replace_all(.data[[group_var]], "_", " ")
      )

    # Base plot
    p <- ggplot2::ggplot(tax_abund, ggplot2::aes(x = .data$Sample, y = .data$Abundance, fill = .data$Taxon)) +
      ggplot2::geom_col(width = arglist$bar_width, color = "black") +
      ggplot2::scale_fill_manual(tax_rank, values = col_vector) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = NULL,
        y = y_lab,
        title = plot_title
      ) +
      theme_microEDA(
        graph_type = "bar",
        sample_text = arglist$sample_text,
        strip_angle = arglist$strip_angle,
        legend_text = arglist$legend_text,
        legend_key = arglist$legend_key,
        show_samples = show_samples
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

    # Apply faceting
    if (is.null(facet_var)) {
      p <- p +
        ggplot2::facet_grid(
          cols = ggplot2::vars(.data[[group_var]]),
          scales = "free_x",
          space = "free",
          switch = "x", # If "x", the top labels will be displayed at the bottom, if NULL at top
        )
    } else {
      p <- p +
        ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]], .data[[group_var]]), scales = "free_x") +
        ggplot2::theme(strip.text.x = ggtext::element_markdown(angle = 0))
    }
    return(p)
  }
}


#' Plot heatmap of taxon abundance or prevalence
#'
#' Creates a ggplot2-based heatmap that displays either mean abundance or
#' prevalence of taxa in its gradient across samples or sample groups.
#' It helps identify dominant or prevalent taxa and visualize patterns possibly
#' influenced by experimental factors.
#' Text labels within cells display mean abundance ± SD and prevalence.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU and taxonomic data.
#' @param tax_rank `Character` string specifying taxonomic rank (e.g., "Genus", "Family").
#'   Abbreviations (e.g., "g") are allowed.
#' @param group_var (Optional) `Character` string indicating a sample variable to group samples.
#'   If `NULL`, all samples are treated as one group.
#' @param plot_title (Optional) `Character` string for the plot title.
#' @param group_labels (Optional) Named character vector to rename group levels (e.g., `c("Old" = "New")`).
#' @param color_by `Character`. What to color the tiles by: `"abundance"` (mean abundance) or `"prevalence"`.
#' @param prevalence_display `Character`. Display prevalence as `"percent"` (e.g., 75%) or `"fraction"` (e.g., 6/8).
#' @param order_by `Character`. Order taxa by `"abundance"` (mean abundance), `"prevalence"`, or `"alphabetical"`.
#' @param min_abundance `Numeric` value. Minimum abundance threshold for a feature to be retained.
#'   Must be non-negative. Features with abundance below this are considered absent.
#' @param min_prevalence `Numeric` value. Minimum prevalence required for retention.
#'   If value is < 1, interpreted as proportion of samples; otherwise, as absolute number of samples.
#' @param as_relative `Logical`. If `TRUE`, applies relative abundance transformation (TSS)
#'   to the input counts. If `FALSE`, uses raw counts. (Default: `TRUE`).
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in plotting without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param ... Additional arguments for fine-tuning. Can include:
#'   \itemize{
#'     \item `k`: `Integer`. Number of decimal places for label formatting (Default: 1).
#'     \item `font_size`: `Numeric` value. Font size for text labels heatmap cells (Default: 2).
#'     \item `legend_text`: Font size for legend (Default: 9).
#'     \item `legend_key`: Size/height of legend (Default: 15).
#'     \item `ntaxa`: Number of top taxa to display (Default: 30).
#'     \item `abundance_criterion`: `Character` string. Criterion to use for filtering:
#'   \describe{
#'     \item{\code{prevalence}:}{Retain features present in at least `min_prevalence` samples
#'       (within group if `group_var` is used) and with abundance >= `min_abundance`
#'       in those samples.}
#'     \item{\code{mean}:}{Also requires that the mean abundance across samples (or group)
#'       is >= `min_abundance`.}
#'   }
#'   Default: `"prevalence"`.
#'     \item `group_requirement`: `Character` string. When `group_var` is specified, determines
#'   whether the filter criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#'     \item `keep_filtered`: Whether to keep filtered out taxa as "Other" (Default: `TRUE`). Takes effect only for `microEDA` objects.
#'     \item `rm_missing`: `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae"). (Default: `FALSE`).
#'     \item `process_taxon`: Whether to clean taxon names (remove underscores, return in italic). (Default: `TRUE`).
#'     \item `low_color`, `high_color`: Colours for low and high ends of the gradient. (Default: `"white"` and `"red"`)
#'   }
#'
#' @details
#' - The function prepares the taxonomic profile by:
#'   \itemize{
#'     \item Filtering taxa by `min_abundance`, `min_prevalence`.
#'     \item Aggregating counts by `tax_rank`.
#'     \item Reducing the total number of taxa to no more than `ntaxa`.
#'     \item Optionally applies relative abundance transformation (TSS) to the data.
#'   }
#' - If `group_var` is specified, samples are grouped by that variable.
#' - Sample counts per group are appended to group labels (e.g., "**Control (n=10)**").
#'
#' The value used for coloring (`color_by`) and the value used for ordering
#' (`order_by`) can be specified independently. Text labels within cells display
#' mean abundance ± SD and prevalence.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' mpa <- microEDA(merged_metaphlan_profiles)
#' plot_taxa_heatmap(mpa, "Species", prevalence_display = "fraction")
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete labs theme element_blank
#' @importFrom dplyr inner_join mutate group_by summarize n n_distinct
#' @importFrom dplyr if_else rename
#' @importFrom tibble tibble rownames_to_column
#' @importFrom rlang sym
#' @importFrom glue glue
#' @importFrom forcats fct_reorder fct_inorder
#' @importFrom stats sd
#'
#' @export
plot_taxa_heatmap <- function(me,
                              tax_rank,
                              group_var = NULL,
                              plot_title = NULL,
                              group_labels = NULL,
                              color_by = c("abundance", "prevalence"),
                              prevalence_display = c("percent", "fraction"),
                              order_by = c("abundance", "prevalence", "alphabetical"),
                              min_abundance = 0,
                              min_prevalence = 0,
                              as_relative = TRUE,
                              filter_by_group = FALSE,
                              ...) {
  # Input validation
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

  if (!is.null(group_var) && !is.character(group_var)) {
    stop("`group_var` must be a character string or NULL.")
  }

  color_by <- match.arg(color_by, choices = c("abundance", "prevalence"))
  order_by <- match.arg(order_by, choices = c("abundance", "prevalence", "alphabetical"))
  prevalence_display <- match.arg(prevalence_display, choices = c("percent", "fraction"))

  if (as_relative) {
    transform <- "TSS"
  } else {
    transform <- "None"
  }

  defaults <- list(
    k = 1,
    font_size = 2,
    legend_text = 9,
    legend_key = 15,
    ntaxa = 30,
    abundance_criterion = "prevalence",
    group_requirement = "any", # Only takes effect when group_var != NULL
    keep_filtered = TRUE,
    rm_missing = FALSE,
    process_taxon = TRUE,
    low_color = "#FFFFFF",
    high_color = "#FF0000"
  )
  # Handle optional ellipsis arguments
  arglist <- list(...)

  arglist <- .warn_invalid_args(allowed = names(defaults), arglist = arglist)

  arglist <- utils::modifyList(defaults, arglist)

  metadata <- NULL
  if (!is.null(group_var)) {
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
  }

  tax_abund <- .prepare_tax_profile(me,
    tax_rank = tax_rank,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    ntaxa = arglist$ntaxa,
    group_var = group_var,
    abundance_criterion = arglist$abundance_criterion,
    filter_by_group = filter_by_group,
    group_requirement = arglist$group_requirement,
    keep_filtered = arglist$keep_filtered,
    rm_missing = arglist$rm_missing,
    transform = transform,
    add_prefix = FALSE,
    process_taxon = arglist$process_taxon,
    calculate_prevalence = TRUE
  )

  is_rel_abund <- "Relative_Abundance" %in% colnames(tax_abund)

  abund_lab <- ifelse(is_rel_abund, "Relative<br>Abundance<br>(%)", "Abundance")

  if (is_rel_abund) {
    tax_abund <- dplyr::rename(tax_abund, "Abundance" = "Relative_Abundance")
  }

  if (is.null(group_var)) {
    n_samples <- dplyr::n_distinct(tax_abund$Sample)

    tax_abund <- tax_abund |>
      dplyr::group_by(.data$Taxon) |>
      dplyr::summarize(
        Samples = paste0("**Samples** (n=", n_samples, ")"),
        mean_abund = mean(.data$Abundance),
        sd = sd(.data$Abundance),
        prevalence = dplyr::first(.data$Prevalence),
        prevalence_n = dplyr::first(.data$Prevalence_n),
        .groups = "drop"
      )
    group_var <- "Samples"
  } else {
    tax_abund <- tax_abund |>
      mutate("Group" = tidyr::replace_na(.data$Group, "Unknown")) |>
      group_by(.data$Group) |>
      mutate(num_samples = dplyr::n_distinct(.data$Sample)) |>
      mutate(
        "Group" = glue::glue("**{.data$Group}** (n={.data$num_samples})"),
        "Group" = str_replace_all(.data$Group, "_", " ")
      )

    tax_abund <- tax_abund |>
      dplyr::group_by(.data$Group, .data$Taxon) |>
      dplyr::summarize(
        num_samples = dplyr::n(),
        mean_abund = mean(.data$Abundance),
        sd = sd(.data$Abundance),
        prevalence = dplyr::first(.data$Prevalence),
        prevalence_n = dplyr::first(.data$Prevalence_n),
        .groups = "drop"
      )
    group_var <- "Group"
  }

  fill_value <- if (color_by == "abundance") "mean_abund" else "prevalence"
  fill_label <- if (color_by == "abundance") paste0("Mean<br>", abund_lab) else "Prevalence<br>(%)"

  # Apply ordering
  if (order_by == "abundance") {
    tax_abund$Taxon <- forcats::fct_reorder(tax_abund$Taxon, tax_abund$mean_abund, .desc = FALSE)
  } else if (order_by == "prevalence") {
    tax_abund$Taxon <- forcats::fct_reorder(tax_abund$Taxon, tax_abund$prevalence, .desc = FALSE)
  } else {
    tax_abund$Taxon <- forcats::fct_rev(tax_abund$Taxon)
  }

  message("Plotting data...")
  p <- ggplot2::ggplot(tax_abund, ggplot2::aes(x = .data[[group_var]], fill = !!rlang::parse_expr(fill_value), y = .data$Taxon)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(
      name = fill_label,
      low = arglist$low_color, high = arglist$high_color,
      expand = c(0, 0),
      limits = c(0, NA),
      # guide = ggplot2::guide_colorbar(barheight = unit(1.5, "cm")) # Adjust legend height
    ) + # add 0 to legend
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      x = NULL, # hide axis labels
      y = NULL,
      title = plot_title
    ) +
    theme_microEDA(
      graph_type = "heat",
      legend_text = arglist$legend_text,
      legend_key = arglist$legend_key
    )

  if (prevalence_display == "percent") {
    p <- p + ggplot2::geom_text(
      aes(label = paste0(
        .specify_decimal(.data$mean_abund, k = arglist$k),
        " \u00B1 ",
        .specify_decimal(sd, k = arglist$k),
        " (", .specify_decimal(.data$prevalence, k = arglist$k), "%)"
      )),
      size = arglist$font_size, vjust = 0.5, hjust = 0.5
    )
  } else {
    p <- p + ggplot2::geom_text(
      aes(label = paste0(
        .specify_decimal(.data$mean_abund, k = arglist$k),
        " \u00B1 ",
        .specify_decimal(sd, k = arglist$k),
        " (", .data$prevalence_n, ")"
      )),
      size = arglist$font_size, vjust = 0.5, hjust = 0.5
    )
  }
  return(p)
}


#' Create an UpSet Plot of Shared and Unique Taxa
#'
#' Generates a ggplot2-based UpSet plot to visualize the intersection patterns of
#' taxa across different sample groups.
#' This provides a scalable alternative to Venn diagrams for microbiome data.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU table, taxonomy, and sample data.
#' @param tax_rank `Character` string specifying the taxonomic rank to plot
#'   (e.g., "Phylum", "Family"). Must be a valid rank in the taxonomy table.
#' @param group_var `Character` string indicating a sample metadata variable
#'   to group samples.
#' @param plot_title Optional `character` string for the plot title. If `NULL`,
#'   a default title is generated based on `tax_rank`.
#' @param show_names `Logical`. If `TRUE`, displays taxon names on intersection
#'   bars with italicized formatting. Default is `FALSE`.
#' @param group_labels Named character `vector` mapping old group names to
#'   new labels (e.g., `c("Old" = "New")`).
#' @param min_abundance `Numeric` value. Minimum abundance threshold for a feature to be retained.
#'   Must be non-negative. Features with abundance below this are considered absent.
#' @param min_prevalence `Numeric` value. Minimum prevalence required for retention.
#'   If value is < 1, interpreted as proportion of samples; otherwise, as absolute number of samples.
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in plotting without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param text_size `Numeric`. Font size for taxon labels when
#'   `show_names = TRUE`. Default is 2.
#' @param ... Additional arguments for fine-tuning. Can include:
#'   \itemize{
#'     \item `abundance_criterion`: `Character` string. Criterion to use for filtering:
#'   \describe{
#'     \item{\code{prevalence}:}{Retain features present in at least `min_prevalence` samples
#'       (within group if `group_var` is used) and with abundance >= `min_abundance`
#'       in those samples.}
#'     \item{\code{mean}:}{Also requires that the mean abundance across samples (or group)
#'       is >= `min_abundance`.}
#'   }
#'   Default: `"prevalence"`.
#'     \item `group_requirement`: `Character` string. When `group_var` is specified, determines
#'   whether the filter criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#'     \item `keep_filtered`: Whether to keep filtered out taxa as "Other"  (Default: `TRUE`). Takes effect only for `microEDA` objects.
#'     \item `rm_missing`: `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae"). (Default: `FALSE`).
#'   }
#'
#' @return A `ggplot` object.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates the input `me` object.
#'   \item Filters taxa by `min_abundance`, `min_prevalence`.
#'   \item Aggregates counts by `tax_rank`.
#'   \item Determines taxon presence/absence per sample group.
#'   \item Constructs an UpSet plot using the `ComplexUpset` package.
#' }
#'
#' Taxon names are formatted in italics when `show_names = TRUE`.
#'
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' plot_taxa_upset(GlobalPatterns, "Phylum", "SampleType")
#'
#' @importFrom rlang sym
#' @importFrom dplyr arrange inner_join mutate filter group_by n_distinct reframe
#' @importFrom tidyr replace_na pivot_wider all_of
#' @importFrom glue glue
#' @importFrom stringr str_replace str_replace_all str_detect str_c
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ylab
#'
#' @export
plot_taxa_upset <- function(me,
                            tax_rank,
                            group_var,
                            plot_title = NULL,
                            show_names = FALSE,
                            group_labels = NULL,
                            min_abundance = 0,
                            min_prevalence = 0,
                            filter_by_group = FALSE,
                            text_size = 2,
                            ...) {
  # Input validation
  if (!inherits(me, "phyloseq")) {
    stop("'me' must be a microEDA or phyloseq object.")
  }

  if (missing(group_var) || is.null(group_var)) {
    stop("'group_var' can't be missing! A variable of interest for the intersections is needed.")
  }

  if (!.is_counts(otu_table(me), silent = TRUE) && !.is_proportion(otu_table(me), silent = TRUE)) {
    stop("'otu_table' is neither counts nor relative abundance - data was likely transformed.")
  }

  if (!.is_valid_rank(tax_rank)) {
    stop(.valid_ranks_msg)
  }
  tax_rank <- .get_full_tax_rank(tax_rank) # In case it was abbreviated

  defaults <- list(
    abundance_criterion = "prevalence",
    group_requirement = "any", # Only takes effect when group_var != NULL
    keep_filtered = TRUE,
    rm_missing = FALSE
  )
  # Handle optional ellipsis arguments
  arglist <- list(...)

  arglist <- .warn_invalid_args(allowed = names(defaults), arglist = arglist)

  arglist <- utils::modifyList(defaults, arglist)

  if (.check_sample_data(me)) {
    if (!is.null(group_var) && !(group_var %in% names(phyloseq::sample_data(me)))) {
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
    abundance_criterion = arglist$abundance_criterion,
    filter_by_group = filter_by_group,
    group_requirement = arglist$group_requirement,
    keep_filtered = arglist$keep_filtered,
    rm_missing = arglist$rm_missing,
    transform = "None",
    add_prefix = FALSE,
    process_taxon = FALSE,
    calculate_prevalence = FALSE
  )

  if (is.null(plot_title) && !is.null(tax_rank)) {
    plot_title <- paste0("Shared and unique ", .get_full_tax_rank(tax_rank, return_plural = TRUE))
  }

  tax_presence <- tax_abund |>
    dplyr::arrange(.data$Sample) |>
    dplyr::inner_join(metadata, by = "Sample") |>
    dplyr::mutate(!!group_var_sym := tidyr::replace_na(.data[[group_var]], "Unknown")) |>
    dplyr::group_by(.data[[group_var]], .data$Sample) |>
    dplyr::filter(.data$Abundance > min_abundance) |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::mutate(num_samples = dplyr::n_distinct(.data$Sample)) |>
    dplyr::mutate(
      !!group_var_sym := glue::glue("{.data[[group_var]]} (n={.data$num_samples})"),
      !!group_var_sym := stringr::str_replace_all(.data[[group_var]], "_", " ")
    ) |>
    dplyr::reframe("Taxon" = unique(.data$Taxon)) |>
    dplyr::group_by(.data$Taxon) |>
    dplyr::mutate(presence = TRUE) |>
    tidyr::pivot_wider(names_from = tidyr::all_of(group_var), values_from = "presence", values_fill = FALSE)

  condition <- setdiff(names(tax_presence), "Taxon")

  if (show_names) {
    italic_label <- tax_presence |>
      dplyr::mutate(
        # Handle '_unclassified' - capture group before it
        Taxon = stringr::str_replace(.data$Taxon, "^(.*)_unclassified$", "Unclass\\.~italic('\\1')"),
        # Taxon = stringr::str_replace(Taxon, "^(.*)_unclassified$", "Unclassified\nitalic('\\1')"),
        # Handle 'Unclassified <name>' -> 'Unclass.~italic(name)'
        Taxon = stringr::str_replace(.data$Taxon, "^Unclassified\\s+(.*)$", "Unclass.~italic('\\1')"),
        # Wrap any remaining plain names (that don't already have ~italic) in italic()
        # to prevent nested italic. This is to ensure "Unclass."~italic('...') stays "flat"
        Taxon = dplyr::if_else(
          !stringr::str_detect(.data$Taxon, "~italic\\("),
          stringr::str_c("italic('", .data$Taxon, "')"),
          .data$Taxon
        ),
        # Replace underscores with ~ for spacing in plotmath
        Taxon = stringr::str_replace_all(.data$Taxon, "_", "~"),
        # Collapse multiple tildes into one
        Taxon = stringr::str_replace_all(.data$Taxon, "~+", "~"),
        # Avoid leading/trailing tildes
        Taxon = stringr::str_replace_all(.data$Taxon, "^~|~$", ""),
        # Special case: revert "italic('Other')" to just 'Other'
        Taxon = stringr::str_replace(
          .data$Taxon,
          "^italic\\('Other([^']*)'\\)$",
          "'Other\\1'"
        )
      )

    ComplexUpset::upset(tax_presence,
      condition,
      name = plot_title,
      set_sizes = ComplexUpset::upset_set_size() + ggplot2::ylab("Total Taxa"),
      base_annotations = list(
        "Intersection size" = (
          ComplexUpset::intersection_size(
            bar_number_threshold = 1,
            color = "grey9",
            fill = "grey80"
          )
          + ggplot2::geom_text(
              mapping = ggplot2::aes(label = italic_label$Taxon),
              position = ggplot2::position_stack(),
              na.rm = TRUE,
              vjust = 1,
              size = text_size,
              parse = TRUE
            )
        )
      ),
      width_ratio = 0.15,
      height_ratio = 0.2
    )
  } else {
    ComplexUpset::upset(tax_presence,
      condition,
      name = plot_title,
      set_sizes = ComplexUpset::upset_set_size() + ggplot2::ylab("Total Taxa"),
      base_annotations = list(
        "Intersection size" = (
          ComplexUpset::intersection_size(
            bar_number_threshold = 1,
            color = "grey9",
            fill = "grey80"
          )
        )
      ),
      width_ratio = 0.15,
      height_ratio = 0.2
    )
  }
}


#' Create a Sankey Plot of Taxonomic Abundances
#'
#' Visualizes the flow of taxa across hierarchical taxonomic ranks (e.g.,
#' Kingdom -> Phylum -> Class -> Order -> Family -> Genus -> Species) using a Sankey diagram.
#' If multiple samples are provided or `samples = NULL` (default), mean abundances
#' across samples are used. The number of displayed taxa is limited to `ntaxa`
#' to avoid visual clutter.
#'
#' @param me A `microEDA` or `phyloseq` object containing OTU table, tax_table, and optionally sample_data.
#' @param tax_rank `Character` string specifying the taxonomic rank to agglomerate at
#'   (e.g., "genus", "family").
#' @param samples Optional `character` vector of sample names to include. If `NULL`, all samples are used.
#' @param group_var Optional name of a variable in `sample_data(me)` to group samples.
#'   Used for filtering and labeling if `filter_by_group = TRUE`.
#' @param plot_title Optional `character` string for the plot title. If `NULL`,
#'   a default title is generated based on the number of samples.
#' @param rm_missing `Logical.` If `TRUE`, removes taxa with missing/unclassified entries
#'   at the specified rank. If `FALSE`, fills missing values by propagating the
#'   last known ancestor, labeling them as "Unclassified _Last_Known_Parent_Clade_"
#'   (e.g., "Unclassified Enterobacteriaceae"). (Default: `FALSE`).
#' @param group_labels Named character `vector` mapping old group names to
#'   new labels (e.g., `c("Old" = "New")`).
#' @param min_abundance `Numeric` value. Minimum abundance threshold for a feature to be retained.
#'   Must be non-negative. Features with abundance below this are considered absent.
#' @param min_prevalence `Numeric` value. Minimum prevalence required for retention.
#'   If value is < 1, interpreted as proportion of samples; otherwise, as absolute number of samples.features.
#' @param as_relative `Logical`. If `TRUE`, applies relative abundance transformation (TSS)
#'   to the input counts. If `FALSE`, uses raw counts. (Default: `TRUE`).
#' @param filter_by_group `Logical`. If `FALSE` (default), filtering is applied
#'   globally across all samples even if `group_var` is specified. This allows
#'   using `group_var` for stratification in plotting without affecting the
#'   filtering scope. If `TRUE`, filtering is applied within each group
#'   defined by `group_var` and `group_requirement`. See [filter_features] for more details on filtering arguments.
#' @param text_size Label size of taxa names (Default: 3).
#' @param ntaxa Number of top taxa to display before grouping into "Other" (Default: 30).
#' @param ... Additional arguments for fine-tuning. Can include:
#'   \itemize{
#'     \item `abundance_criterion`: `Character` string. Criterion to use for filtering:
#'   \describe{
#'     \item{\code{prevalence}:}{Retain features present in at least `min_prevalence` samples
#'       (within group if `group_var` is used) and with abundance >= `min_abundance`
#'       in those samples.}
#'     \item{\code{mean}:}{Also requires that the mean abundance across samples (or group)
#'       is >= `min_abundance`.}
#'   }
#'   Default: `"prevalence"`.
#'     \item `group_requirement`: `Character` string. When `group_var` is specified, determines
#'   whether the filter criterion must be met in `"any"` group or `"all"` groups.
#'   Default: `"any"`.
#'     \item `keep_filtered`: Whether to keep filtered out taxa as "Other"  (Default: `TRUE`). Takes effect only for `microEDA` objects.
#'     \item `ncp`: Number of control points on the Bezier curve that forms the edge.
#'     Larger numbers will result in smoother curves, but cost more computational time. (Default: 100).
#'     \item `slope`: Slope parameter for the Bezier curves used to depict the edges.(Default: 0.5).
#'   }
#'
#' @return A `ggplot` object displaying the Sankey diagram.
#'
#' @examples
#' data(GlobalPatterns, package = "phyloseq")
#' plot_taxa_sankey(GlobalPatterns, tax_rank = "Order", samples = "Even3")
#'
#' @importFrom ggsankeyfier geom_sankeynode geom_sankeyedge pivot_stages_longer StatSankeynode
#' @importFrom ggplot2 ggplot aes labs scale_fill_viridis_c scale_x_discrete
#' @importFrom ggtext geom_richtext
#' @importFrom stringi stri_replace_all_regex stri_replace_first_regex
#' @importFrom tidyr pivot_longer
#'
#' @export
plot_taxa_sankey <- function(me,
                             tax_rank,
                             samples = NULL,
                             group_var = NULL,
                             plot_title = NULL,
                             rm_missing = FALSE,
                             group_labels = NULL,
                             min_abundance = 0,
                             min_prevalence = 0,
                             as_relative = TRUE,
                             filter_by_group = FALSE,
                             text_size = 3,
                             ntaxa = 30,
                             ...) {
  # Input validation
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

  if (!tax_rank %in% phyloseq::rank_names(me)) {
    stop(
      "Rank '", tax_rank, "' not present in tax_table. Available ranks: ",
      paste(phyloseq::rank_names(tax_rank), collapse = ", ")
    )
  }
  rank_index <- which(phyloseq::rank_names(me) == tax_rank)
  tax_ranks <- phyloseq::rank_names(me)[1:rank_index]

  if (as_relative) {
    transform <- "TSS"
    abundance_type <- "Relative\nAbundance\n(%)"
  } else {
    transform <- "None"
    abundance_type <- "Abundance"
  }

  defaults <- list(
    abundance_criterion = "prevalence",
    group_requirement = "any", # Only takes effect when group_var != NULL
    keep_filtered = TRUE,
    slope = 0.5,
    ncp = 100
  )
  # Handle optional ellipsis arguments
  arglist <- list(...)

  arglist <- .warn_invalid_args(allowed = names(defaults), arglist = arglist)

  arglist <- utils::modifyList(defaults, arglist)

  arglist$abundance_criterion <- match.arg(arglist$abundance_criterion, choices = c("prevalence", "mean"))
  arglist$group_requirement <- match.arg(arglist$group_requirement, choices = c("any", "all"))

  if (!is.null(group_var)) {
    if (.check_sample_data(me)) {
      if (!(group_var %in% names(phyloseq::sample_data(me)))) {
        stop("'group_var' not found in sample metadata.")
      }

      if (!is.null(group_labels)) {
        phyloseq::sample_data(me) <- .rename_values(phyloseq::sample_data(me), group_var, group_labels)
      }
    } else {
      stop("Can't use 'group_var' without sample metadata!")
    }
  }

  group_var_arg <- if (filter_by_group) group_var else NULL

  if (!is.null(samples)) {
    samples <- intersect(phyloseq::sample_names(me), samples)
    if (length(samples) < 1) stop("No matching sample(s) found.")

    me <- phyloseq::prune_samples(samples, me)

    if (is.null(plot_title)) {
      if (length(samples) > 1) {
        plot_title <- paste0("Mean Relative Abundance for samples:", " *", toString(samples), "*")
      } else {
        plot_title <- paste0("Relative Abundance for sample:", " *", samples, "*")
      }
    }
  } else if (is.null(plot_title)) {
    plot_title <- "Mean Relative Abundance over all samples"
  }

  # Temporarily enable immediate warning print so it appears before agglomeration warnings
  old_warn <- options(warn = 1)
  on.exit(options(old_warn), add = TRUE)

  me <- filter_features(me,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    group_var = group_var_arg,
    abundance_criterion = arglist$abundance_criterion,
    group_requirement = arglist$group_requirement,
    keep_filtered = arglist$keep_filtered
  )

  me <- agglomerate_taxa(me,
    tax_rank = tax_rank,
    rm_missing = rm_missing,
    transform = transform,
    add_prefix = FALSE
  )

  rank_index <- which(phyloseq::rank_names(me) == tax_rank)
  tax_ranks <- phyloseq::rank_names(me)[1:rank_index]

  # Calculate the mean per Taxon if there is more than 1 sample
  if (phyloseq::nsamples(me) > 1) {
    abundance_type <- paste0("Mean\n", abundance_type)

    merged <- data.frame(tax_table(me), phyloseq::otu_table(me))

    sample_cols <- names(which(vapply(merged, is.numeric, logical(1L))))

    merged_long <- pivot_longer(merged, cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance")

    mean_by_rank <- lapply(tax_ranks, function(rank_col) {
      merged_long |>
        # CRITICAL STEP: Sum / aggregate abundances of all OTUs belonging to the
        # same Taxon WITHIN each Sample, otherwise the mean will reflect the
        # mean abundance within OTUs not the Taxon
        dplyr::group_by(dplyr::across(dplyr::all_of(rank_col)), .data$Sample) |>
        dplyr::summarise(Sample_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") |>
        # Now calculate the mean of summed abundances across samples
        dplyr::group_by(dplyr::across(dplyr::all_of(rank_col))) |>
        dplyr::summarise(mean_abund = mean(.data$Sample_Abundance, na.rm = TRUE), .groups = "drop") |>
        dplyr::rename(Taxon = dplyr::all_of(rank_col)) |>
        dplyr::mutate(Rank = rank_col)
    })

    mean_all_ranks <- do.call(rbind, mean_by_rank)


    Abundance <- as.data.frame(phyloseq::otu_table(me))
    Abundance <- rowMeans(Abundance)
    merged <- data.frame(tax_table(me), Abundance)

    # Add abundance to Taxon label

    # 1. Convert mean_all_ranks to a named vector for fast lookup by Taxon and Rank
    mean_lookup <- mean_all_ranks |>
      dplyr::select(.data$Taxon, .data$mean_abund) |>
      tibble::deframe() |>
      split(mean_all_ranks$Rank) # Nested list: Rank -> Taxon -> mean_abund

    # # 2. Function to append mean abundance to a taxon if found
    # append_mean <- function(x, rank) {
    #   if (is.na(x) || !is.character(x)) {
    #     return(x)
    #   }
    #   mean_val <- mean_lookup[[rank]][[x]]
    #   if (!is.null(mean_val)) {
    #     return(paste0(x, " (", round(mean_val, 2), ")"))
    #   } else {
    #     return(x)
    #   }
    # }

    # 3. Apply across relevant columns with known ranks
    for (col in names(merged)) {
      if (col %in% tax_ranks) {
        # rank_name <- col
        # merged[[col]] <- vapply(merged[[col]], append_mean, character(1L), rank = rank_name, USE.NAMES = FALSE)

        # Replace the vapply loop with a vectorized approach using dplyr::case_when
        # Extract mean values for this rank
        rank_means <- mean_all_ranks |>
          dplyr::filter(.data$Rank == col) |>
          dplyr::select(.data$Taxon, .data$mean_abund) |>
          tibble::deframe()

        # Append mean abundance to taxon labels
        merged <- merged |>
          dplyr::mutate(
            !!sym(col) := {
              current_val <- .data[[col]]
              dplyr::case_when(
                # Check if value exists, is not "Other", and exists in rank_means
                !is.na(current_val) & current_val != "Other" & current_val %in% names(rank_means) ~
                  paste0(current_val, " (", round(rank_means[current_val], 2), ")"),
                # Otherwise keep value as is
                TRUE ~ current_val
              )
            }
          )
      }
    }
  } else { # For single sample
    Abundance <- as.data.frame(phyloseq::otu_table(me)) |>
      dplyr::rename(Abundance = all_of(samples))
    merged <- data.frame(tax_table(me), Abundance)

    # Add abundance to Taxon label
    for (col in tax_ranks) {
      # 1. Calculate sums per group
      sums <- merged |>
        dplyr::group_by(dplyr::across(dplyr::all_of(col))) |>
        dplyr::summarise(temp_sum = sum(.data$Abundance, na.rm = TRUE), .groups = "drop")

      # 2. Join back to original data
      merged <- merged |>
        dplyr::left_join(sums, by = col) |>
        # 3. Update the column with "Name (Sum)"
        # mutate(across(all_of(col), ~ paste0(.x, " (", round(temp_sum, 2), ")"))) |>
        dplyr::mutate(dplyr::across(dplyr::all_of(col), ~ ifelse(grepl("^Other *\\(<.*%\\)$", .x),
          .x, # Keep original
          paste0(.x, " (", round(.data$temp_sum, 2), ")")
        ))) |>
        dplyr::select(-.data$temp_sum)
    }
  }

  merged_reduced <- .apply_incrementing_filter(merged, ntaxa = ntaxa, tax_column = tax_rank, increment = 0.1)
  merged <- merged_reduced$data

  other_idx <- which(rownames(merged) == "Other")
  # Replace NA values in the "Other" row with "Other (<...%) from tax_rank
  merged[other_idx, ][is.na(merged[other_idx, ])] <- merged[other_idx, tax_rank]

  merged <- merged |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.character),
        ~ .x |>
          # 0. Clean slate: Remove ANY existing asterisks to prevent doubling
          # stringi::stri_replace_all_regex("\\*", "") |>
          # 1. Remove rank prefixes (k__, p__, etc.)
          stringi::stri_replace_all_regex("[kpcofgst]__", "") |>
          # 2. Replace underscores with spaces
          stringi::stri_replace_all_regex("_", " ") |>
          # 3. Handle 'unclassified' WITH a name (e.g., "Unclassified Firmicutes (5.2)")
          # (?i): Case-insensitive flag
          stringi::stri_replace_first_regex(
            "(?i)^unclassified\\s+([A-Za-z].*)\\s+\\(([^)]+)\\)$",
            "**Unclass.** ***$1*** ($2)"
          ) |>
          # 4. Handle generic terms like "Unclassified" or "Other" WITHOUT a specific name
          # These should become **Unclassified** (...) or **Other** (<...%)
          stringi::stri_replace_first_regex(
            "(?i)^(unclassified|Other|Unknown)\\s+\\(([^)]+)\\)$",
            "**$1** ($2)"
          ) |>
          # 5. Standard formatting for specific named taxa (e.g., "Actinomycetia (17.37)")
          # Use negative lookahead for any double asterisk, so replacement is only
          # applied if string does NOT start with "**" (already formatted)
          stringi::stri_replace_first_regex(
            "^(?!\\*\\*)(.*)\\s+\\(([^)]+)\\)$",
            "***$1*** ($2)"
          )
      )
    )

  merged_long <- ggsankeyfier::pivot_stages_longer(merged, stages_from = tax_ranks, values_from = "Abundance")

  pos <- ggsankeyfier::position_sankey(order = "ascending", align = "center", v_space = "auto")

  ggplot2::ggplot(merged_long, ggplot2::aes(
    x = .data$stage, y = Abundance, group = .data$node,
    connector = .data$connector, edge_id = .data$edge_id
  )) +
    ggsankeyfier::geom_sankeynode(position = pos) +
    ggsankeyfier::geom_sankeyedge(
      position = pos, aes(fill = Abundance),
      colour = "black",
      linewidth = 1,
      ncp = arglist$ncp,
      slope = arglist$slope
    ) +
    ggtext::geom_richtext(
      aes(label = .data$node),
      size = text_size,
      # stat = "sankeynode",
      stat = ggsankeyfier::StatSankeynode,
      position = pos,
      fill = NA, label.color = NA,
      label.padding = unit(rep(0, 4), "pt")
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = plot_title
    ) +
    ggplot2::scale_fill_viridis_c(abundance_type, option = "turbo") +
    scale_x_discrete(expand = ggplot2::expansion(add = c(0.3, .5))) +
    theme_microEDA(graph_type = "sankey")
}
