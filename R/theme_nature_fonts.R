#' Custom Minimal Font Theme for Nature-style Plots
#'
#' A lightweight ggplot2 theme that standardizes font sizes for all common 
#' plot elements. Ideal for consistent presentation in publication-style figures.
#'
#' @param base_size Numeric. Base font size used for all text elements.
#'
#' @return A ggplot2 theme object.
#' @export
theme_nature_fonts <- function(base_size = 5) {
  ggplot2::theme(
    axis.text     = ggplot2::element_text(size = base_size),
    axis.title    = ggplot2::element_text(size = base_size),
    plot.title    = ggplot2::element_text(size = base_size, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5),
    legend.title  = ggplot2::element_text(size = base_size),
    legend.text   = ggplot2::element_text(size = base_size),
    strip.text    = ggplot2::element_text(size = base_size),
    plot.caption  = ggplot2::element_text(size = base_size, hjust = 0)
  )
}