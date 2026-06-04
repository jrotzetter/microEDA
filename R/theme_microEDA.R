#' Custom Theme Function for Microbiome Data Plots
#'
#' A modular ggplot2 theme function tailored for microbiome composition plots,
#' with specific styling for bar, heatmap, and sankey plots.
#'
#' @param graph_type `Character` string specifying the plot type.
#'   One of `NULL`, `"bar"`, `"heat"`, or `"sankey"`.
#'   Controls which theme elements are modified.
#' @param sample_text Font size for x-axis sample labels (used in bar plots).
#' @param strip_angle Angle for facet strip text labels (bar plots only).
#' @param legend_text Font size for legend text.
#' @param legend_key Size of legend keys in points.
#' @param show_samples `Logical`. If `FALSE`, hides x-axis sample labels (bar plots only).
#'
#' @return A `ggplot2::theme` object built upon `theme_classic()`.
#'
#' @details
#' - **Bar plots** (`graph_type = "bar"`):
#'   - Rotates x-axis text 90 degrees.
#'   - Adds horizontal y-grid lines (`gray60`).
#'   - Centers plot title.
#'   - Formats legend with markdown support via `ggtext`.
#'   - Optionally hides sample labels.
#' - **Heatmaps** (`graph_type = "heat"`):
#'   - Enables markdown in axis and legend text.
#'   - Removes axis lines and ticks for cleaner layout.
#' - **Sankey plots** (`graph_type = "sankey"`):
#'   - Hides y-axis elements (title, text, line, ticks).
#'   - Applies markdown to plot title.
#'
#' @importFrom ggplot2 theme_classic theme element_text element_blank unit
#' @importFrom ggtext element_markdown
#'
#' @keywords internal
#' @noRd
theme_microEDA <- function(graph_type = NULL,
                           sample_text = 8,
                           strip_angle = 0,
                           legend_text = 10,
                           legend_key = 7,
                           show_samples = TRUE) {
  base_theme <- ggplot2::theme_classic()

  if (graph_type == "bar") {
    microEDA_theme <- base_theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggtext::element_markdown(size = legend_text),
        legend.key.size = ggplot2::unit(legend_key, "pt"),
        axis.text.x = ggplot2::element_text(angle = 90, size = sample_text, vjust = 0.5),
        strip.placement = "outside",
        strip.text.x = ggtext::element_markdown(angle = strip_angle),
        panel.grid.major.y = ggplot2::element_line(colour = "gray60", linewidth = 0.8)
      )
    if (!show_samples) {
      microEDA_theme <- microEDA_theme + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
      )
    }
  }

  if (graph_type == "heat") {
    microEDA_theme <- base_theme +
      ggplot2::theme(
        axis.text.x.top = ggtext::element_markdown(vjust = 0.5),
        axis.text.y = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown(size = 8),
        axis.line = ggplot2::element_blank(), # hide axis line
        axis.ticks = ggplot2::element_blank(), # hide axis ticks
        legend.text = ggplot2::element_text(size = legend_text),
        legend.key.height = ggplot2::unit(legend_key, "pt")
      )
  }

  if (graph_type == "sankey") {
    microEDA_theme <- base_theme +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.title = ggtext::element_markdown()
      )
  }
  return(microEDA_theme)
}
