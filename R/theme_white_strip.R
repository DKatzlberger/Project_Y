#' White Strip Theme
#'
#' A ggplot2 theme component that sets a white background for facet strips,
#' places them inside the plotting area, and uses black strip text.
#'
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme object.
#' @export
theme_white_strip <- function(...) {
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.text       = ggplot2::element_text(color = "black"),
    strip.placement  = "inside",
    ...
  )
}
