#' Plot Sample Proportions with Labels per Group and Class
#'
#' @param data A data.frame with grouping columns.
#' @param x Character. Column name for x-axis (e.g., ancestry).
#' @param fill Character. Column name for fill (e.g., class/subtype).
#' @param point_size Numeric. Controls the size of the text labels and base font size in the plot.
#'
#' @return A ggplot object.
#' @export
plot_output_column_proportion <- function(
  data,
  x,
  fill,
  point_size = 0.5
) {
  ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x    = !!rlang::sym(x),
      fill = !!rlang::sym(fill)
      )
   ) +
    ggplot2::geom_bar(
      position = "fill"
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = scales::percent(
          ggplot2::after_stat(count) / sum(ggplot2::after_stat(count)),
          accuracy = 0.1
        )
      ),
      stat     = "count",
      position = ggplot2::position_fill(vjust = 0.5),
      color    = "black",
      size     = point_size * 3
    ) +
    ggplot2::labs(
      y = "Proportion"
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background() +
    theme_small_legend(
      base_size = point_size * 10
    ) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}
