#' Plot Correlation Bar Plots with Facets and Annotations
#'
#' @param data A data.frame with columns: method, train_test, train_infer, difference, p_param
#' @param test_label Character. Label for test column (default = "Test").
#' @param infer_label Character. Label for inference column (default = "Inference").
#' @param point_size Numeric. Controls annotation text and base font size (default = 0.5).
#'
#' @return A ggplot2 object.
#' @export
plot_correlation <- function(
  data,
  test_label   = "Test",
  infer_label  = "Inference",
  point_size   = 0.5,
  title        = NULL
) {

  # Reshape to long format for ggplot
  long_data <- reshape2::melt(
    data,
    id.vars       = "method",
    measure.vars  = c("train_test", "train_infer"),
    variable.name = "phase",
    value.name    = "correlation"
  )

  # Relabel phases for x-axis
  long_data$phase <- factor(
    long_data$phase,
    levels = c("train_test", "train_infer"),
    labels = c(test_label, infer_label)
  )

  # Merge annotations per method
  annotation_data <- unique(data[, c("method", "difference", "p_param")])

  # Create the plot
  ggplot2::ggplot(
    data    = long_data,
    mapping = ggplot2::aes(
      x     = phase,
      y     = correlation
    )
  ) +
    ggplot2::geom_col(
      width    = point_size
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(method)
    ) +
    ggplot2::geom_text(
      data    = annotation_data,
      mapping = ggplot2::aes(
        x     = -Inf,
        y     = Inf,
        label = paste0(
          "Delta = ", signif(difference, 3), ", ",
          "p = ", signif(p_param, 3)
        )
      ),
      hjust        = -0.05,
      vjust        = 1.1,
      size         = point_size * 3,
      inherit.aes  = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x     = NULL,
      y     = "Correlation"
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
      legend.position = "none"
    )
}
