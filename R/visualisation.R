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
#'     \item `rm_missing`: Whether to remove taxa with missing names  (Default: `FALSE`).
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
#' plot_composition(mpa, "Species")
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
plot_composition <- function(me,
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
  stopifnot(inherits(me, "phyloseq"))

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
#'     \item `rm_missing`: Whether to remove taxa with missing names (Default: `FALSE`).
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
#' plot_heatmap(mpa, "Species", prevalence_display = "fraction")
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
plot_heatmap <- function(me,
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
  stopifnot(inherits(me, "phyloseq"))

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
