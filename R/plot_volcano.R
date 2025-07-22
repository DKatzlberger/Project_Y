plot_volcano <- function(
  data,
  effect_col,         # Column with effect size (e.g., logFC)
  p_value_col,        # Column with raw p-values
  effect_label,       # Label for x-axis (e.g., "log₂ Fold Change")
  sig_thr    = 0.05,  # Significance threshold (p-value)
  effect_thr = 1,     # Effect size threshold
  point_size = 0.5,   # Size of points
  y_cap      = 5,     # Max value for -log10(p)
  title      = NULL,  # Plot title
  color_by   = "significant"  # "significant" or a column like "prop_signif"
) {
  # Compute log p-values and significance flags
  data$neg_log_p <- -log10(data[[p_value_col]])
  data$significant <- abs(data[[effect_col]]) > effect_thr & data[[p_value_col]] < sig_thr

  # Cap the y-values
  data$neg_log_p_capped <- pmin(data$neg_log_p, y_cap)
  data$y_capped <- data$neg_log_p > y_cap

  # Base plot setup
  p <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x     = !!rlang::sym(effect_col),
      y     = neg_log_p_capped,
      color = !!rlang::sym(color_by)
    )
  ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    # Add triangle markers for capped points
    ggplot2::geom_point(
      data = data[data$y_capped, ],
      aes(
        x = !!rlang::sym(effect_col),
        y = y_cap
      ),
      shape = 25,
      fill  = "white",
      size  = point_size
    ) +
    # Significance threshold lines
    ggplot2::geom_hline(
      yintercept = -log10(sig_thr),
      color      = "blue",
      linetype   = "dashed",
      linewidth  = point_size / 1.5
    ) +
    ggplot2::geom_vline(
      xintercept = c(-effect_thr, effect_thr),
      color      = "blue",
      linetype   = "dashed",
      linewidth  = point_size / 1.5
    ) +
    ggplot2::labs(
      title = title,
      x     = effect_label,
      y     = paste0("-log10(", p_value_col, ")")
    )

  # Color logic
  if (color_by == "significant") {
    p <- p +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "red", "FALSE" = "gray"),
        name   = "Significant"
      ) +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- p +
      ggplot2::scale_color_viridis_c(
        name   = "Reproducibility",
        option = "D"
      ) +
      ggplot2::theme(legend.position = "right")
  }

  p <- p +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_small_legend(
      base_size = point_size * 10
    ) +
    theme_white_background()

  return(p)
}
