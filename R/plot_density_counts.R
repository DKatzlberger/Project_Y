plot_density_of_counts <- function(
  data_matrix,
  n_samples = 12,
  x_axis_label,
  point_size = 0.5
) {
  # Subset the first n samples (rows)
  data_subset <- data_matrix[1:n_samples, , drop = FALSE]

  # Add sample (patient) labels
  data_df <- as.data.frame(data_subset)
  data_df$Patient <- rownames(data_subset)

  # Reshape to long format using base R
  data_long <- reshape(
    data_df,
    varying = colnames(data_subset),
    v.names = "Value",
    timevar = "Gene",
    times = colnames(data_subset),
    idvar = "Patient",
    direction = "long"
  )

  # Plot using ggplot2
  p <- ggplot(
    data_long,
    aes(
      x = Value,
      color = Patient
    )
  ) +
    geom_density(
      fill = NA
    ) +
    labs(
      x = x_axis_label,
      y = "Density"
    ) +
    guides(
      color = guide_legend(ncol = 3)
    ) +
    theme_nature_fonts(base_size = point_size * 10) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend(base_size = point_size * 10) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )

  return(p)
}
