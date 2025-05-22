#' Plot Mean-Variance Trend (Voom-like)
#'
#' Generates a mean-variance trend plot to visualize the relationship between
#' mean gene expression and variability across samples. Uses a square root
#' transformation of standard deviation to stabilize variance.
#'
#' Assumes data is in samples x features (e.g., genes) format.
#'
#' @param data A numeric matrix or data frame with samples in rows and genes in columns.
#'             Expression values should be non-negative.
#' @param x_axis_label Character. Label for the x-axis (typically the type of transformed expression).
#' @param point_size Numeric. Size of points in the plot. Default is 0.5.
#'
#' @return A ggplot object showing the mean vs. square-root standard deviation per gene.
#' @export
plot_mean_variance_trend <- function(
  data,
  x_axis_label,
  point_size = 0.5
) {

  # Compute per-gene mean and SD across samples
  gene_means <- colMeans(data, na.rm = TRUE)
  gene_sds   <- apply(data, 2, sd, na.rm = TRUE)

  df <- data.frame(
    mean_raw    = gene_means,
    sqrt_sd_raw = sqrt(gene_sds)
  )

  ggplot2::ggplot(
    data    = df, 
    mapping = ggplot2::aes(
        x = mean_raw, 
        y = sqrt_sd_raw
      )
    ) +
    ggplot2::geom_point(
        size = point_size
    ) +
    ggplot2::geom_smooth(
      formula   = y ~ x,
      method    = "loess",
      color     = "red",
      se        = FALSE,
      linewidth = point_size 
    ) +
    ggplot2::labs(
      x = x_axis_label,
      y = "Sqrt(Standard Deviation)"
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background()
}
