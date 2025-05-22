#' Plot a Volcano Plot with Significance and Effect Thresholds
#'
#' @param data A data.frame containing statistical test results.
#' @param effect_col Character. Name of the effect size column (e.g., log2 fold change).
#' @param p_value_col Character. Name of the p-value column.
#' @param effect_label Character. X-axis label (default = "Effect Size").
#' @param significance_thresh Numeric. P-value significance threshold (default = 0.05).
#' @param effect_thresh Numeric. Absolute effect size threshold (default = 1).
#' @param point_alpha Numeric. Point transparency (default = 0.6).
#' @param point_size Numeric. Controls the size of the points (default = 1).
#' @param title Character. Plot title.
#'
#' @return A ggplot2 object.
#' @export
plot_volcano <- function(
  data,
  effect_col,
  p_value_col,
  effect_label,
  sig_thr    = 0.05,
  effect_thr = 1,
  point_size = 0.5,
  title      = NULL
) {

  data$neg_log_p   <- -log10(data[[p_value_col]])
  data$significant <- abs(data[[effect_col]]) > effect_thr & data[[p_value_col]] < sig_thr

  ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x     = !!rlang::sym(effect_col),
      y     = neg_log_p,
      color = significant
    )
  ) +
    ggplot2::geom_point(
      size  = point_size
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "gray"),
      name   = "Significant"
    ) +
    ggplot2::labs(
      title = title,
      x     = effect_label,
      y     = paste0("-log10(", p_value_col, ")")
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background() +
    theme_small_legend(
      base_size = point_size * 10
    ) +
    ggplot2::theme(
      legend.position = "none"
    )
}
