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
  return(microEDA_theme)
}
