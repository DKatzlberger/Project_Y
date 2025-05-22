# THEMES
# Font theme
theme_nature_fonts <- function(base_size = 5) {
  theme(
    axis.text     = element_text(size = base_size),
    axis.title    = element_text(size = base_size),
    plot.title    = element_text(size = base_size, hjust = 0.5),
    plot.subtitle = element_text(size = base_size, hjust = 0.5),
    legend.title  = element_text(size = base_size),
    legend.text   = element_text(size = base_size),
    strip.text    = element_text(size = base_size),
    plot.caption  = element_text(size = base_size, hjust = 0)
  )
}

# Small legend theme
theme_small_legend <- function(...) {
  theme(
    legend.key.spacing = unit(0, "cm"),
    legend.key.height  = unit(0.3, "cm"),  
    legend.key.width   = unit(0.3, "cm"),
    legend.margin      = margin(0, 0, 0, 0),
    ...
  )
}

# White background
theme_white_background <- function(...) {
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    panel.grid.major = element_blank(),                           
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.5, "mm"),                    
    ...
  )
}

# White strip background
theme_white_strip <- function(...) {
  theme(
    strip.background = element_rect(fill = "white", color = NA),  
    strip.text = element_text(color = "black"),   
    strip.placement = "inside",                                  
    ...
  )
}

# Zero margin
theme_zero_margin <- function(...){
  theme(
    plot.margin = margin(0, 0, 0, 0),  
    ...
  )
}

# PLOTS
# Save function
save_ggplot <- function(plot, save_path, width = 5, height = 5, with_legend = TRUE) {
  # Try to save the plot
  tryCatch({
    # Check if the plot is NULL or not a ggplot object
    if (is.null(plot) || !inherits(plot, "gg")) {
      warning("The plot variable is not assigned or is invalid. The plot will not be saved.")
      return(NULL)  # Exit the function gracefully
    }

    # Optionally remove the legend if requested
    if (!with_legend) {
      plot <- plot + theme(legend.position = "none")  # Remove legend if needed
    }

    # Save the plot to a file
    ggsave(filename = save_path, plot = plot, width = width, height = height)
    
    # Remove the plot from the global environment using the name of the variable
    plot_name <- deparse(substitute(plot))
    rm(list = plot_name, envir = .GlobalEnv)

  }, error = function(e) {
    # If an error occurs, catch it and show a message without stopping the script
    warning("Error in creating/saving the plot: ", conditionMessage(e))
    return(NULL)  # Exit gracefully on error
  })
  
  invisible()  # Ensures nothing is returned after deletion
}

# QC
plot_output_column_proportion <- function(data, x, fill) {
  # Plot
  p <- ggplot(
    data = data,
    aes(
        x    = !!sym(x), 
        fill = !!sym(fill)
      )
  ) +
  geom_bar(
    position = "fill"
  ) +
  geom_text(
    aes(
      label = scales::percent(after_stat(count)/sum(after_stat(count)), accuracy = 0.1)
    ),
    stat = "count",
    position = position_fill(vjust = 0.5),  
    color = "black",  
    size = 2
  ) +
  labs(
    y = "Proportion",
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend() +
  theme(
    legend.position = "bottom"
  )
  
  return(p)
}

plot_output_column_count <- function(data, x, fill) {
  # Plot
  p <- ggplot(
    data = data,
    aes(
        x    = !!sym(x), 
        fill = !!sym(fill)
      )
  ) +
  geom_bar(
    position = "dodge"
  ) +
  geom_text(
    aes(
      label = after_stat(count)
    ),
    stat = "count", 
    position = position_dodge(width = 0.8), 
    vjust = -0.5,  
    color = "black", 
    size = 2     
  ) +
  labs(
    y    = "Count",
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend() +
  theme(
    legend.position = "bottom"
  )
  
  return(p)
}

plot_density_of_samples <- function(data_matrix, n_samples = 12, x_axis_label) {

  # Extract the first n samples (rows) from the matrix
  data_subset <- data_matrix[1:n_samples, ]
  
  # Convert the matrix to a data frame and reshape it to long format using native pipe
  data_long <- data_subset |>
    as.data.frame() |>
    rownames_to_column("Patient") |>
    pivot_longer(
      cols      = -Patient, 
      names_to  = "Gene", 
      values_to = "Value"
    )
  
  # Plot density for the first n samples
  p <- ggplot(
    data_long, 
    aes(
      x     = Value, 
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
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
  
  return(p)
}

plot_qq_of_genes <- function(data_matrix, n_features = 10) {

  # Convert matrix to a long format (samples x genes)
  expression_df <- as.data.frame(data_matrix) |>
    pivot_longer(
      cols      = everything(), 
      names_to  = "Gene", 
      values_to = "Expression"
    ) 

  # Filter features
  expression_df <- expression_df |>
      filter(
        Gene %in% names(as.data.frame(data_matrix))[1:n_features]
        )

  # Calculate Pearson correlation for each gene (Sample Quantiles vs Theoretical Quantiles)
  correlation_values <- expression_df |>
    group_by(Gene) |>
    summarize(
      correlation = cor(Expression, qqnorm(Expression, plot.it = FALSE)$x),
      .groups     = 'drop'
    )
  
  # Merge
  expression_df <- expression_df |>
    left_join(
      correlation_values, 
      by = "Gene"
    )

  p <- ggplot(
    expression_df, 
    aes(
      sample = Expression
    )
  ) +
    stat_qq(
      size  = 0.5,
      shape = 1
    ) +
    stat_qq_line(
      linewidth = (0.5 / 2)
    ) +
    geom_text(
      aes(
        x     = -Inf, 
        y     = Inf, 
        label = paste("Pearson: ", round(correlation, 3)),
        hjust = -0.1, 
        vjust = 1.1
      ), 
      size     = 2,
      fontface = "plain"
    ) +
    facet_grid(
      cols = vars(Gene)
    ) +
    labs(
      x = "Theoretical Quantiles", 
      y = "Sample Quantiles"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip()
  
  return(p)
}

# Mean variance trend plot (Voom like plot)
plot_mean_variance_trend <- function(data, x_axis_label) {
  # Generate a mean-variance trend plot (Voom-like plot).
  
  # This function visualizes the relationship between the mean and variance of gene expression
  # levels across samples. It assumes that **samples are rows** and **genes are columns**.

  # The function computes the mean and standard deviation of each gene across samples,
  # applies a square root transformation to the standard deviation, and plots the trend
  # using a LOESS smoother.

  # Args:
  #     data (matrix or data frame): A numeric matrix or data frame where **rows represent samples**
  #                                  and **columns represent genes**. Expression values should be non-negative.
  #     x_axis_label (str): Label for the x-axis, typically describing the mean expression per gene.

  # Returns:
  #     ggplot object: A mean-variance trend plot showing the relationship between mean expression and variance.

  
  # Step 1: Compute mean and standard deviation across samples (columns)
  mean_raw <- colMeans(data)  # Mean for each gene
  sd_raw   <- apply(data, 2, sd)  # Standard deviation for each gene
  
  # Step 2: Create data frame for plotting
  df <- data.frame(
    mean_raw    = mean_raw,
    sqrt_sd_raw = sqrt(sd_raw)  
  )

  # Step 3: Create the plot
  p <- ggplot(
      df, 
      aes(
        x = mean_raw, 
        y = sqrt_sd_raw
      )
    ) +
    geom_point(size = 0.5) +  # Raw data points
    geom_smooth(
      formula   = y ~ x,
      method    = "loess", 
      color     = "red", 
      se        = FALSE, 
      linewidth = 0.5
    ) +  
    labs(
      x = x_axis_label,  
      y = "Sqrt(Standard Deviation)"
    ) +
    theme_nature_fonts() +
    theme_white_background() 
  
  return(p)
}

# DESCRIPTIVE MODEL BUILDING

plot_variable_count <- function(data, var, point_size = 0.5) {

  # Add facet column
  data <- mutate(
    data,
    facet = "Count"
  )

  # Count unique values in the variable
  n_levels <- n_distinct(data[[var]])

  # Base plot
  p <- ggplot(
    data, 
    aes(
      x    = !!sym(var), 
      fill = !!sym(var)
    )
  ) +
    geom_bar(position = "dodge") +
    facet_grid(cols = vars(facet)) +
    labs(
      title = var,
      x     = NULL,
      y     = NULL,
      fill  = var
    ) +
    theme_nature_fonts(base_size = (point_size * 10)) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1)
    )

  # Add legend columns if more than 8 levels
  if (n_levels > 8) {
    p <- p + guides(fill = guide_legend(ncol = 2))
  }

  return(p)
}

plot_variable_proportion <- function(data, var, point_size = 0.5) {

  # Add facet column and dummy x
  data <- mutate(
    data,
    facet   = "Proportion",
    dummy_X = var
  )

  # Count unique values in the variable
  n_levels <- n_distinct(data[[var]])

  # Base plot
  p <- ggplot(
    data, 
    aes(
      x    = dummy_X, 
      fill = !!sym(var)
    )
  ) +
    geom_bar(position = "fill") +
    facet_grid(cols = vars(facet)) +
    labs(
      x    = NULL,
      y    = NULL, 
      fill = var
    ) +
    theme_nature_fonts(base_size = (point_size * 10)) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1)
    )

  # Add legend columns if more than 8 levels
  if (n_levels > 8) {
    p <- p + guides(fill = guide_legend(ncol = 2))
  }

  return(p)
}

plot_variable_hist <- function(data, var, point_size = 0.5){

  # Convert the column to numeric and overwrite it
  data <- data |>
    mutate(
      !!sym(var) := as.numeric(.data[[var]]),
      facet       = "Count"
    ) |>
    filter(!is.na(.data[[var]]))

  # Plot
  p <- ggplot(
    data = data,
    aes(
      x = !!sym(var)
      )
    ) +
    geom_histogram(
      bins = 30
    ) +
    facet_grid(
      cols = vars(facet)
    ) +
    labs(
      title = var,
      x     = NULL,
      y     = NULL
    ) + 
    theme_nature_fonts(
      base_size = (point_size * 10)
    ) +
    theme_white_background() +
    theme_white_strip()

  return(p)
}

plot_clusters <- function(data, x, y, color, title, point_size = 0.5){

  # Plot
  p <- ggplot(
    data,
    aes(
      x     = !!sym(x),
      y     = !!sym(y),
      color = !!sym(color)
    )
  ) +
  geom_point(
    size = point_size
  ) +
  labs(
    title = title
  ) +
  theme_nature_fonts(
    base_size = (point_size * 10)
  ) +
  theme_white_background() +
  theme_small_legend()

  return(p)
}

plot_pc_associated <- function(data, x, point_size = 0.5){

  p <- ggplot(
  data = data,
  aes(
    x = !!sym(x),
    y = r2
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = sig_feature
    ),
    vjust = -0.1,  
    size = 3
  ) +
  facet_wrap(
    ~ variable,
    ncol = 3  
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) + 
  labs(
    x = "Principal Component",
    y = "Variance explained (R2)"
  ) +
  theme_nature_fonts(
    base_size = (point_size * 10)
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() 

  return(p)
}

plot_pc_variance <- function(data, point_size = 0.5){

  p  <- ggplot(
  data = data, 
  aes(
    x = pc, 
    y = pc_r2
    )
  ) + 
  geom_col() +
  labs(
    x = "Principal Component",
    y = "Variance explained (R2)"
  ) +
  theme_nature_fonts(
    base_size = (point_size * 10)
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() 

  return(p)
}

plot_pc_combinations <- function(data, color_by, pc_pairs, point_size = 0.5) {
  # Create a map from PC to label
  pc_label_map <- data %>%
    select(pc, pc_label) %>%
    distinct() %>%
    deframe()

  # Pivot to wide format
  coords_wide <- data %>%
    select(idx, pc, coordinates) %>%
    pivot_wider(names_from = pc, values_from = coordinates)

  # Add metadata
  meta <- data %>%
    select(idx, all_of(color_by)) %>%
    distinct()

  coords_wide <- left_join(coords_wide, meta, by = "idx")

  # Create individual plots
  plots <- map(pc_pairs, function(pair) {
    pc_x <- pair[1]
    pc_y <- pair[2]

    label_x <- pc_label_map[pc_x]
    label_y <- pc_label_map[pc_y]

    df <- coords_wide %>%
      select(x = all_of(pc_x), y = all_of(pc_y), all_of(color_by))

    ggplot(
      df, 
      aes(
        x = x, 
        y = y, 
        color = .data[[color_by]]
        )
      ) +
      geom_point(
        size = point_size
      ) +
      labs(
        title = color_by,
        x     = label_x,
        y     = label_y,
        color = color_by
      ) +
      theme_nature_fonts(
        base_size = (point_size * 10)
      ) +
      theme_white_background() +
      theme_white_strip() +
      theme_small_legend()
    }
  )

  # Combine plots
  combined_plot <- wrap_plots(plots, ncol = 2) + plot_layout(guides = "collect") & theme(legend.position = "right")

  return(combined_plot)
}



# CROSS ANCESTRY DGE
# Null distributions
plot_cor_diff_histogram <- function(data, x, point_size = 0.5) {
  # Ensure x_var is a symbol (unquoted column name)
  x_sym <- rlang::ensym(x)

  data <- data |>
  mutate(
    Legend = dplyr::case_when(
      Statistic == "H0_bootstrapped" ~ "Bootstrapped differences",
      Statistic == "Observed"        ~ "Observed difference"
    )
  )
  
  p <- ggplot(
    data = data,
    aes(x = !!x_sym)
  ) +
  geom_histogram(
    data = filter(data, Statistic == "H0_bootstrapped"),
    bins = 50
  ) +
  geom_segment(
    data = filter(data, Statistic == "Observed"),
    aes(
      x     = !!x_sym,
      xend  = !!x_sym,
      y     = 0,
      yend  = Inf,
      color = Legend  
    )
  ) +
  facet_grid(
    cols = vars(Metric)
  ) +
  scale_color_manual(
    name   = "Statistic",
    values = c(
      "Observed difference"      = "orange"
    )
  ) +
  labs(
    x     = "Difference of correlation coefficient to train set (Subset - Ancestry)",
    y     = "Count"
  ) +
  theme_nature_fonts(
    base_size = point_size * 10
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  coord_cartesian(clip = "off") 
  
  # Extract histogram data
  built     <- ggplot_build(p)
  hist_data <- built$data[[1]]
  
  # Compute range and offsets
  x_range <- range(hist_data$x, na.rm = TRUE)
  x_min   <- x_range[1]
  x_max   <- x_range[2]
  x_span  <- x_max - x_min
  
  arrow_start_neg <- 0 - 0.01 * x_span
  arrow_end_neg   <- x_min + 0.05 * x_span
  arrow_start_pos <- 0 + 0.01 * x_span
  arrow_end_pos   <- x_max - 0.05 * x_span
  
  max_count <- max(hist_data$count, na.rm = TRUE)
  arrow_y_below <- -max_count * 0.05
  text_y_below  <- arrow_y_below - max_count * 0.05

  # Add arrows and annotations
  p <- p +
    annotate(
      "segment",
      x          = arrow_start_neg, xend = arrow_end_neg,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_neg + arrow_end_neg) / 2,
      y     = text_y_below,
      label = "Ancestry closer",
      size  = point_size * 3.5,
      hjust = 0.5
    ) +
    annotate(
      "segment",
      x          = arrow_start_pos, xend = arrow_end_pos,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_pos + arrow_end_pos) / 2,
      y     = text_y_below,
      label = "Subset closer",
      size  = point_size * 3.5,
      hjust = 0.5
    ) 
  
  return(p)
}

plot_euclid_diff_histogram <- function(data, x, point_size = 0.5) {
  # Ensure x_var is a symbol (unquoted column name)
  x_sym <- rlang::ensym(x)

  data <- data |>
  mutate(
    Legend = dplyr::case_when(
      Statistic == "H0_bootstrapped" ~ "Bootstrapped differences",
      Statistic == "Observed"        ~ "Observed difference"
    )
  )
  
  p <- ggplot(
    data = data,
    aes(x = !!x_sym)
  ) +
  geom_histogram(
    data = filter(data, Statistic == "H0_bootstrapped"),
    bins = 50
  ) +
  geom_segment(
    data = filter(data, Statistic == "Observed"),
    aes(
      x     = !!x_sym,
      xend  = !!x_sym,
      y     = 0,
      yend  = Inf,
      color = Legend  
    )
  ) +
  facet_wrap(
    ~ Feature
  ) +
  scale_color_manual(
    name   = "Statistic",
    values = c(
      "Observed difference" = "orange"
    )
  ) +
  labs(
    x     = "Difference of euclidean distance to train set (Subset - Ancestry)",
    y     = "Count"
  ) +
  theme_nature_fonts(
    base_size = point_size * 10
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  coord_cartesian(clip = "off")
  
  # Extract histogram data
  built     <- ggplot_build(p)
  hist_data <- built$data[[1]]
  
  # Compute range and offsets
  x_range <- range(hist_data$x, na.rm = TRUE)
  x_min   <- x_range[1]
  x_max   <- x_range[2]
  x_span  <- x_max - x_min
  
  arrow_start_neg <- 0 - 0.01 * x_span
  arrow_end_neg   <- x_min + 0.05 * x_span
  arrow_start_pos <- 0 + 0.01 * x_span
  arrow_end_pos   <- x_max - 0.05 * x_span
  
  max_count <- max(hist_data$count, na.rm = TRUE)
  arrow_y_below <- -max_count * 0.05
  text_y_below  <- arrow_y_below - max_count * 0.05

  # Add arrows and annotations
  p <- p +
    annotate(
      "segment",
      x          = arrow_start_neg, xend = arrow_end_neg,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_neg + arrow_end_neg) / 2,
      y     = text_y_below,
      label = "Subset closer",
      size  = point_size * 3.5,
      hjust = 0.5
    ) +
    annotate(
      "segment",
      x          = arrow_start_pos, xend = arrow_end_pos,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_pos + arrow_end_pos) / 2,
      y     = text_y_below,
      label = "Ancestry closer",
      size  = point_size * 3.5,
      hjust = 0.5
    )
  
  return(p)
}

plot_logFC_diff_histogram <- function(data, x, point_size = 0.5) {
  # Ensure x_var is a symbol (unquoted column name)
  x_sym <- rlang::ensym(x)

  data <- data |>
  mutate(
    Legend = dplyr::case_when(
      Statistic == "H0_bootstrapped" ~ "Bootstrapped differences",
      Statistic == "Observed"        ~ "Observed difference"
    )
  )
  
  p <- ggplot(
    data = data,
    aes(x = !!x_sym)
  ) +
  geom_histogram(
    data = filter(data, Statistic == "H0_bootstrapped"),
    bins = 50
  ) +
  geom_segment(
    data = filter(data, Statistic == "Observed"),
    aes(
      x     = !!x_sym,
      xend  = !!x_sym,
      y     = 0,
      yend  = Inf,
      color = Legend  
    )
  ) +
  facet_wrap(
    ~ Feature
  ) +
  scale_color_manual(
    name   = "Statistic",
    values = c(
      "Observed difference" = "orange"
    )
  ) +
  labs(
    x     = "Difference of logFC (Subset - Ancestry)",
    y     = "Count"
  ) +
  theme_nature_fonts(
    base_size = point_size * 10
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  coord_cartesian(clip = "off")
  
  # Extract histogram data
  built     <- ggplot_build(p)
  hist_data <- built$data[[1]]
  
  # Compute range and offsets
  x_range <- range(hist_data$x, na.rm = TRUE)
  x_min   <- x_range[1]
  x_max   <- x_range[2]
  x_span  <- x_max - x_min
  
  arrow_start_neg <- 0 - 0.01 * x_span
  arrow_end_neg   <- x_min + 0.05 * x_span
  arrow_start_pos <- 0 + 0.01 * x_span
  arrow_end_pos   <- x_max - 0.05 * x_span
  
  max_count <- max(hist_data$count, na.rm = TRUE)
  arrow_y_below <- -max_count * 0.05
  text_y_below  <- arrow_y_below - max_count * 0.05

  # Add arrows and annotations
  p <- p +
    annotate(
      "segment",
      x          = arrow_start_neg, xend = arrow_end_neg,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_neg + arrow_end_neg) / 2,
      y     = text_y_below,
      label = "Ancestry upreg.",
      size  = point_size * 3.5,
      hjust = 0.5
    ) +
    annotate(
      "segment",
      x          = arrow_start_pos, xend = arrow_end_pos,
      y          = arrow_y_below, yend = arrow_y_below,
      arrow      = arrow(length = unit(0.15, "cm")),
      linewidth  = 0.3
    ) +
    annotate(
      "text",
      x     = (arrow_start_pos + arrow_end_pos) / 2,
      y     = text_y_below,
      label = "Subset upreg.",
      size  = point_size * 3.5,
      hjust = 0.5
    )
  
  return(p)
}

# Volcano plots
plot_cross_ancestry_volcano <- function(
  data,
  x,
  y,
  x_label,
  arrow_left,
  arrow_right,
  log_y      = TRUE,
  sig        = 0.05,
  thr        = 1,
  point_size = 0.5,
  caption    = NULL
) {
  # Extract values
  x_vals  <- data[[x]]
  y_vals  <- if (log_y) -log10(data[[y]]) else data[[y]]
  y_label <- if (log_y) paste0("-log10(", y,")") else y

  # Symbol
  x_sym <- sym(x)
  y_sym <- sym(y)

  # Main ggplot object
  p <- ggplot(
      data = data,
      aes(
        x     = !!x_sym,
        y     = if (log_y) -log10(!!y_sym) else !!y_sym,
        color = (.data[[y]] < sig & abs(.data[[x]]) > thr)
      )
    ) +
    geom_point(
      size = point_size
    ) +
    scale_color_manual(
      values = c(
        "TRUE" = "red", 
        "FALSE" = "lightgrey"
      )
    ) +
    labs(
      x       = x_label,
      y       = y_label,
      caption = caption
    ) +
    theme_nature_fonts(base_size = point_size * 10) +
    theme_white_background() +
    theme_white_strip() +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off")

  # Compute span for arrow/text positioning
  x_range <- range(x_vals, na.rm = TRUE)
  x_span  <- diff(x_range)
  x_min   <- x_range[1]
  x_max   <- x_range[2]

  arrow_start_neg <- 0 - 0.01 * x_span
  arrow_end_neg   <- x_min + 0.05 * x_span
  arrow_start_pos <- 0 + 0.01 * x_span
  arrow_end_pos   <- x_max - 0.05 * x_span

  y_min <- min(y_vals, na.rm = TRUE)
  arrow_y_below <- y_min - 0.05
  text_y_below  <- arrow_y_below - 0.05

  # Add directional arrows and labels
  p <- p +
    annotate(
      "segment",
      x = arrow_start_neg, xend = arrow_end_neg,
      y = arrow_y_below,   yend = arrow_y_below,
      arrow = arrow(length = unit(0.15, "cm")),
      linewidth = 0.3
    ) +
    # Left error
    annotate(
      "text",
      x     = (arrow_start_neg + arrow_end_neg) / 2,
      y     = text_y_below,
      label = arrow_left,
      size  = point_size * 3.5,
      hjust = 0.5
    ) +
    annotate(
      "segment",
      x = arrow_start_pos, xend = arrow_end_pos,
      y = arrow_y_below,   yend = arrow_y_below,
      arrow = arrow(length = unit(0.15, "cm")),
      linewidth = 0.3
    ) +
    # Arrow right
    annotate(
      "text",
      x     = (arrow_start_pos + arrow_end_pos) / 2,
      y     = text_y_below,
      label = arrow_right,
      size  = point_size * 3.5,
      hjust = 0.5
    )

  # Calculate y limits for segments (spanning points' range)
  y_min <- min(y_vals, na.rm = TRUE)
  y_max <- max(y_vals, na.rm = TRUE)

  p <- p +
    annotate(
      "segment",
      x = -thr, xend = -thr,
      y = y_min, yend = Inf,
      color     = "blue",
      linetype  = "dashed",
      linewidth = point_size / 2
    ) +
    annotate(
      "segment",
      x = thr, xend = thr,
      y = y_min, yend = Inf,
      color     = "blue",
      linetype  = "dashed",
      linewidth = point_size / 2
    ) +
    geom_hline(
      yintercept = if (log_y) -log10(sig) else sig,
      linetype   = "dashed",
      color      = "blue",
      linewidth  = point_size / 2
    ) 

  return(p)
}

plot_boxplot <- function(
  data,
  x,
  y,
  fill,
  x_label,
  point_size = 0.5
) {
  #
  x_sym    <- sym(x)
  y_sym    <- sym(y)
  fill_sym <- sym(fill)

  p <- ggplot(
    data = data,
    aes(
      x    = !!x_sym,
      y    = !!y_sym,
      fill = !!fill_sym
    )
  ) +
  geom_boxplot(
    width        = point_size,
    outlier.size = point_size / 10,
  ) +
  facet_wrap(
    ~ Feature,  
    scales         = "free",  
    strip.position = "top",
    ncol           = 5,
    nrow           = 2
  ) +
  labs(
    x    = x_label,
    y    = "z-score"
  ) +
  theme_nature_fonts(
    base_size = (point_size * 10)
  ) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank()
  )

  return(p)
}


# INTERACTIONS
# Volcano plot
volcano_plot <- function(
  data, 
  x            = "logFC", 
  y            = "adj.P.Val", 
  logFC_thr    = 1, 
  point_size   = 0.5, 
  alpha_value  = 0.5,
  top_features = NULL, 
  facet_rows   = NULL, 
  facet_cols   = NULL
  ) {

  # Ensure top_features is handled correctly
  if (is.null(top_features)) {
    top_features <- character(0) 
  }
  
  # Convert character vectors to quosures for multiple facet variables
  row_facet <- if (!is.null(facet_rows)) rlang::syms(facet_rows) else NULL
  col_facet <- if (!is.null(facet_cols)) rlang::syms(facet_cols) else NULL

  x_sym <- sym(x)
  y_sym <- sym(y)

  p <- data |>
    ggplot(
      aes(
        x     = !!x_sym,
        y     = -log10(!!y_sym),
        color = (!!y_sym < 0.05 & abs(!!x_sym) > logFC_thr)
      )
    ) +
    geom_point(
      size = point_size
    ) +
    geom_text_repel(
      data = data |> filter(Feature %in% top_features),
      aes(
        label = Feature
        ),
      size               = 1.5,
      color              = "black",  
      segment.color      = "black",  
      min.segment.length = 0
    ) +
    facet_grid(
      rows = if (!is.null(row_facet)) vars(!!!row_facet) else NULL, 
      cols = if (!is.null(col_facet)) vars(!!!col_facet) else NULL
    ) + 
    geom_vline(
      xintercept = c(-logFC_thr, logFC_thr), 
      linetype   = "dashed", 
      color      = "blue",
      linewidth  = (point_size/2),
      alpha      = alpha_value
    ) +
    geom_hline(
      yintercept = -log10(0.05), 
      linetype   = "dashed", 
      color      = "blue",
      linewidth  = (point_size/2),
      alpha      = alpha_value
    ) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "lightgrey")
    ) +
    theme_nature_fonts(
      base_size = (point_size * 10)
    ) +
    theme_white_background() +
    theme_white_strip() +
    theme(
      legend.position = "none"
    )

  return(p)
}

ma_plot <- function(data, logFC_thr = 1, y_axis_label, point_size = 0.5, alpha_value = 0.5,
                    top_features = NULL, facet_rows = NULL, facet_cols = NULL) {

  if (is.null(top_features)) {
    top_features <- character(0)
  }

  # Facet quosures
  row_facet <- if (!is.null(facet_rows)) rlang::syms(facet_rows) else NULL
  col_facet <- if (!is.null(facet_cols)) rlang::syms(facet_cols) else NULL

  # Determine significance
  data <- data |>
    dplyr::mutate(significant = adj.P.Val < 0.05 & abs(logFC) > logFC_thr)

  p <- ggplot(
    data, 
    aes(
      x = AveExpr, 
      y = logFC
      )
    ) +
    # Plot non-significant first (background)
    geom_point(
      data  = dplyr::filter(data, !significant),
      color = "lightgrey",
      size  = point_size,
    ) +
    # Then plot significant points on top (foreground)
    geom_point(
      data  = dplyr::filter(data, significant),
      color = "red",
      size  = point_size,
    ) +
    # Annotate top features
    geom_text_repel(
      data = dplyr::filter(data, Feature %in% top_features),
      aes(
        label = Feature
        ),
      size               = 1.5,
      color              = "black",
      segment.color      = "black",
      min.segment.length = 0
    ) +
    # Faceting
    facet_grid(
      rows = if (!is.null(row_facet)) vars(!!!row_facet) else NULL,
      cols = if (!is.null(col_facet)) vars(!!!col_facet) else NULL
    ) +
    # Threshold lines
    geom_hline(
      yintercept = c(-logFC_thr, logFC_thr),
      linetype   = "dashed",
      color      = "blue",
      linewidth  = (point_size / 2),
      alpha      = alpha_value
    ) +
    # Axis labels and styling
    labs(
      x = paste("Mean", y_axis_label)
    ) +
    theme_nature_fonts(
      base_size = (point_size * 10)
    ) +
    theme_white_background() +
    theme_white_strip() +
    theme(legend.position = "none")

  return(p)
}



# Venn diagram
venn_diagram <- function(venn_data, title = NULL, base_size) {

  region_data <- venn_data$regionData

  venn_data_intersections <- region_data |>
    filter(count > 0 & grepl("/", name))

  # Adjust coordinates of labels
  venn_data_adjusted <- venn_setlabel(venn_data) |>
    mutate(
      adjusted_X = case_when(
        X == max(X) ~ X * 0.8,  
        X == min(X) ~ X * 0.8,  
        TRUE ~ X  
      ),
      adjusted_Y = case_when(
        X == max(X) ~ Y * 0.8,  
        X == min(X) ~ Y * 1.2,  
        TRUE ~ Y 
      )
    )
  
  # Create the plot and store it in variable 'p'
  p <- ggplot() +
    geom_polygon(
      data = venn_regionedge(venn_data_intersections),  # Use only intersection data
      aes(X, Y, group = id, fill = count),  # Fill by count of the intersection
      color = "black",  # Set a border color to distinguish the regions
      size = 0.5
    ) +
    # 1. region count layer
    # geom_polygon(
    #   data = venn_regionedge(venn_data),
    #   aes(X, Y, group = id, fill = (1/count)),
    #   # fill = "white"
    # ) +
    # 2. set edge layer
    geom_path(
      data = venn_setedge(venn_data), 
      aes(X, Y, group = id), 
      show.legend = FALSE
    ) +
    # 3. set label layer
    geom_text(
      data = venn_data_adjusted,
      aes(adjusted_X, Y, label = name),
      size = 2
    ) +
    # 4. region label layer
    geom_label(
      data = venn_regionlabel(venn_data),
      aes(X , Y, label = count), 
      size = 2
    ) +
    coord_equal() +
    theme_void() +
    scale_fill_gradient(low = "white", high = "blue") + 
    theme(legend.position = "none")
  
  # If title is provided, add it
  if (!is.null(title)) {
      p <- p + ggtitle(title) + 
        theme(
          plot.title = element_text(size = base_size, hjust = 0.5)  
        )
    }
  
  # Return the plot object
  return(p)
}



# Interaction boxplots
interactions_boxplot <- function(expression, x, fill, point_size = 0.5) {
  
  # Generate the plot
  p <- expression |>
    ggplot(
      aes(
        x    = fct_rev(toupper(!!sym(x))),
        y    = zscore,
        fill = !!sym(fill)
      )
    ) +
    geom_boxplot(
      width        = 0.5,
      outlier.size = 0.1,
    ) +
    scale_y_continuous(
      breaks = scales::breaks_extended(n = 3)
    ) +
    facet_wrap(
      ~ Feature,  
      scales         = "free",  
      strip.position = "top",
      ncol           = 5,
      nrow           = 2
    ) +
    labs(
      x = "Ancestry",
      y = "z-score",
      fill = "Condition"
    ) +
    theme_nature_fonts(
      base_size = (point_size * 10)
    ) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank()
    )
  
  # Return the plot object
  return(p)
}

interaction_heatmap <- function(expression, output_column, ancestry_column, path_to_save_location) {
  
  # Expression matrix
  expression_matrix <- expression |>
    dplyr::select(Feature, idx, zscore) |>
    tidyr::spread(key = idx, value = zscore) |>
    tibble::column_to_rownames("Feature") |>
    as.matrix()

  # Ancestry
  train_ancestry <- levels(expression[[ancestry_column]])[1]
  infer_ancestry <- levels(expression[[ancestry_column]])[2]

  
  # Condition
  class_0    <- levels(expression[[output_column]])[1]
  class_1    <- levels(expression[[output_column]])[2]
  comparison <- c(class_0, class_1)
  # Colors
  class_0_color <- "#027c58"
  class_1_color <- "purple"
  colors        <- c(class_0_color, class_1_color)
  # Vector
  condition_colors <- setNames(colors, comparison)
  
  # Annotation
  annotations <- expression |>
    dplyr::select(idx, dplyr::all_of(c(output_column, ancestry_column))) |>
    dplyr::distinct() |>
    dplyr::arrange(match(idx, colnames(expression))) |>
    dplyr::mutate(
      group = paste(.data[[ancestry_column]], .data[[output_column]], sep = "."),
      group = factor(
        group, 
        levels = c(
          paste(train_ancestry, class_0, sep = "."),  
          paste(train_ancestry, class_1, sep = "."),  
          paste(infer_ancestry, class_0, sep = "."),    
          paste(infer_ancestry, class_1, sep = ".")     
        )
      )
    )
  
  # Heatmap annotation
  heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Ancestry = ComplexHeatmap::anno_empty(
      border = FALSE, 
      height = grid::unit(0.2, "cm"), 
      show_name = TRUE
    ),
    Condition = annotations[[output_column]],
    col = list(
      Condition = condition_colors
    ),
    annotation_name_gp = grid::gpar(fontsize = 5),
    border = FALSE,
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    simple_anno_size = grid::unit(0.2, "cm"),
    show_legend = FALSE
  )
  
  # Legends
  ancetsry_lgd <- ComplexHeatmap::Legend(
    title = "Ancestry", 
    at = c(toupper(train_ancestry), toupper(infer_ancestry)), 
    legend_gp = grid::gpar(
      fill = c("#3399ff", "#ff9900")
    ),
    labels_gp = grid::gpar(fontsize = 5),
    title_gp = grid::gpar(fontsize = 5, fontface = "plain"),
    grid_height = grid::unit(0.3, "cm"), 
    grid_width = grid::unit(0.3, "cm"),
    direction = "horizontal"
  )
  
  condition_lgd <- ComplexHeatmap::Legend(
    title = "Condition", 
    at = comparison, 
    legend_gp = grid::gpar(
      fill = c(condition_colors[comparison[1]], condition_colors[comparison[2]])
    ),
    labels_gp = grid::gpar(fontsize = 5),
    title_gp = grid::gpar(fontsize = 5, fontface = "plain"),
    grid_height = grid::unit(0.3, "cm"), 
    grid_width = grid::unit(0.3, "cm"),
    direction = "horizontal"
  )
  
  lgd_list <- list(ancetsry_lgd, condition_lgd)
  
  # Heatmap
  heatmap <- ComplexHeatmap::Heatmap(
    expression_matrix,
    name = "z-score", 
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = annotations$group,
    cluster_column_slices = FALSE,
    column_title = NULL,
    top_annotation = heatmap_annotation,
    row_names_gp = grid::gpar(fontsize = 5),
    row_title_gp = grid::gpar(fontsize = 5),
    row_names_side = "left",
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 5, fontface = "plain"), 
      labels_gp = grid::gpar(fontsize = 5, fontface = "plain"),
      grid_height = grid::unit(0.3, "cm"), 
      grid_width = grid::unit(0.3, "cm"),
      direction = "horizontal"
    )
  )
  
  # Save heatmap
  save_name <- file.path(path_to_save_location, "Interaction_heatmap.pdf")
  grDevices::pdf(save_name, width = 6 , height = 3)
  
  ComplexHeatmap::draw(
    heatmap, 
    annotation_legend_list = lgd_list,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
  
  # Internal helper to draw colored block around ancestry groups
  group_block_anno <- function(group, empty_anno, gp = grid::gpar(), 
                               label = NULL, label_gp = grid::gpar()) {
    seekViewport(glue::glue("annotation_{empty_anno}_{min(group)}"))
    loc1 <- grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
    
    seekViewport(glue::glue("annotation_{empty_anno}_{max(group)}"))
    loc2 <- grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))
    
    seekViewport("global")
    grid::grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
                    just = c("left", "bottom"), gp = gp)
    
    if (!is.null(label)) {
      grid::grid.text(label, x = (loc1$x + loc2$x) * 0.5, y = (loc1$y + loc2$y) * 0.5, gp = label_gp)
    }
  }
  
  # Draw ancestry group background blocks
  group_block_anno(
    1:2, "Ancestry", 
    gp = grid::gpar(fill = "#3399ff", col = NA)
  )
  group_block_anno(
    3:4, "Ancestry", 
    gp = grid::gpar(fill = "#ff9900", col = NA)
  )
  
  invisible(grDevices::dev.off())
}




# Enrichment plot
fgsea_plot <- function(enrichment, x, top_n = NULL, facet_rows = NULL, facet_cols = NULL) {
  
  # Convert character vectors to quosures for multiple facet variables
  row_facet <- if (!is.null(facet_rows)) rlang::syms(facet_rows) else NULL
  col_facet <- if (!is.null(facet_cols)) rlang::syms(facet_cols) else NULL

  # Ensure x_axis is treated as a symbol dynamically
  x_var <- rlang::sym(x)

  # Optionally filter to top_n pathways
  if (!is.null(top_n)) {
    enrichment <- enrichment |> 
      arrange(desc(abs(NES))) |> 
      filter(row_number() <= top_n)
  }

  # Create the plot
  p <- enrichment |>
    ggplot(
      aes(
        x = !!x_var,  
        y = pathway,
        color = NES,
        size = pmin(-log10(padj), 5)  # Size scaled by -log10 of adj.P.Val
      )
    ) +
    geom_point() +
    facet_grid(
      rows = if (!is.null(row_facet)) vars(!!!row_facet) else NULL, 
      cols = if (!is.null(col_facet)) vars(!!!col_facet) else NULL
    ) +
    scale_size_binned(
       range = c(1, 3)    
    ) +
    scale_color_gradient2(
        high = "red", 
        mid = "white", 
        low = "blue"
    ) +
    labs(
      x = x,  # Automatically labels the x-axis
      y = "MSigDB Hallmark 2020 gene set",
      size = "-log10(adj.P.Val)"
    ) +
    theme_nature_fonts() +
    theme_small_legend() +
    theme_white_background() +
    theme_white_strip() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title.position = "top"
    )
  
  # Return the plot
  return(p)
}
