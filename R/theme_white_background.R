#' White Background Theme
#'
#' A minimal ggplot2 theme with a white background, no grid lines, and black axes.
#' Designed for clean, uncluttered visual presentation.
#'
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme object.
#' @export
theme_white_background <- function(...) {
  ggplot2::theme(
    panel.background    = ggplot2::element_rect(fill = "white", color = NA),
    plot.background     = ggplot2::element_rect(fill = "white", color = NA),
    panel.grid.major    = ggplot2::element_blank(),
    panel.grid.minor    = ggplot2::element_blank(),
    axis.line           = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks          = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length   = ggplot2::unit(0.5, "mm"),
    ...
  )
}