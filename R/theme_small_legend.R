#' Compact Legend Theme with Visual Match to Font Size
#'
#' Scales legend keys to visually match the font size, not just by raw pt.
#'
#' @param base_size Numeric. Font size in pt; used for both text and legend scaling.
#' @param ... Additional theme elements.
#'
#' @return A ggplot2 theme.
#' @export
theme_small_legend <- function(base_size = 5, ...) {
  key_size_pt <- base_size * 1.3  # visually balances with text height

  ggplot2::theme(
    legend.key.spacing = ggplot2::unit(0.1, "pt"),
    legend.key.height  = ggplot2::unit(key_size_pt, "pt"),
    legend.key.width   = ggplot2::unit(key_size_pt, "pt"),
    legend.margin      = ggplot2::margin(0, 0, 0, 0),
    ...
  )
}