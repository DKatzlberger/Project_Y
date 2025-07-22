plot_pvalue_concordance <- function(
  data,              # data.table containing p-values from two methods
  p_col_x,           # unadjusted p-value column for method X
  p_col_y,           # unadjusted p-value column for method Y
  signif_col,        # factor column indicating significance source (e.g., "Both", "Limma only", etc.)
  method_x_label,    # label for method X (e.g., "limma")
  method_y_label,    # label for method Y (e.g., "bootstrap")
  log_cap     = 5,   # numeric: maximum -log10(p) value to display
  point_size  = 0.5, # numeric: point size in plot
  title       = NULL
) {

  # Copy
  data <- data.table::copy(data)  

  # Mathematical transformation: logP = -log10(p)
  data[, `:=`(
    logp_x  = -log10(get(p_col_x)),
    logp_y  = -log10(get(p_col_y))
  )]

  # Cap extreme log p-values for visual clarity
  data[, `:=`(
    logp_x_capped = pmin(logp_x, log_cap),
    logp_y_capped = pmin(logp_y, log_cap),
    is_capped_x   = logp_x > log_cap,
    is_capped_y   = logp_y > log_cap,
    capped        = factor(ifelse(logp_x > log_cap | logp_y > log_cap, "Capped", "Uncapped"))
  )]

  # Plot: logP_x vs logP_y with capped visual range
  p <- ggplot(
    data, 
    aes(
      x = logp_x_capped, 
      y = logp_y_capped, 
      color = !!sym(signif_col),
      shape = capped
    )
  ) +
    geom_point(
      size = point_size
    ) +
    geom_abline(
      slope     = 1, 
      intercept = 0, 
      linetype  = "dashed", 
      color     = "blue"
    ) +
    scale_color_manual(
      values = c(
        "None"           = "gray70", 
        "Limma only"     = "#1f77b4",
        "Bootstrap only" = "#ff7f0e",
        "Both"           = "#d62728"
      ),
      name = "Significances source"
    ) +
    scale_shape_manual(
      values = c("Uncapped" = 1, "Capped" = 25),
      name   = "Point Status"
    ) +
    labs(
      title = title,
      x     = paste("-log10(", p_col_x,")"),
      y     = paste("-log10(", p_col_y,")")
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background() +
    theme_small_legend(
      base_size = point_size * 10
    ) +
    theme(legend.position = "right")

  return(p)
}
