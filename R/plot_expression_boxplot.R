#' Plot Z-scored Expression Boxplots with Fixed Color Scheme
#'
#' Generates boxplots of z-scored expression values grouped by a specified variable and colored by another,
#' with one facet per feature. Fill colors are fixed and mapped by the levels of the `fill` column.
#'
#' @param data A data.table or data.frame with z-scored expression and metadata.
#' @param x Character. Column name for the x-axis grouping variable (e.g., ancestry).
#' @param y Character. Column name for the expression value (e.g., z-score).
#' @param fill Character. Column name for the fill color grouping (e.g., subtype or condition).
#' @param feature_column Character. Column name identifying the feature (for faceting). Default is \"feature\".
#' @param point_size Numeric. Controls the size of text and stroke. Default is 0.5.
#'
#' @return A ggplot object.
#' @export
plot_expression_boxplots <- function(
  data,
  x,
  y,
  fill,
  feature_column = "feature",
  point_size     = 0.5
) {
  data.table::setDT(data)

  # Get factor levels of fill column and assign hardcoded colors
  fill_levels <- unique(as.character(data[[fill]]))
  fill_colors <- setNames(c("#027c58", "purple")[seq_along(fill_levels)], fill_levels)

  ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x    = !!rlang::sym(x),
      y    = !!rlang::sym(y),
      fill = !!rlang::sym(fill)
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      linewidth     = point_size / 1.5
    ) +
    ggplot2::geom_jitter(
      mapping = ggplot2::aes(color = !!rlang::sym(fill)),
      shape   = 21,
      size    = point_size,
      alpha   = point_size / 2,
      stroke  = 0,
      position = ggplot2::position_jitterdodge()
    ) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_color_manual(values = fill_colors) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", feature_column)),
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = NULL,
      y = y
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend(
      base_size = point_size * 10
    ) +
    ggplot2::theme(
      legend.position = "right"
    )
}
