# COLORS
genetic_ancestry_colors <- c(
  "admix"= "#ff4d4d", 
  "afr" = "#ff9900", 
  "amr" =  "#33cc33",
  "eur" = "#3399ff", 
  "eas" = "#cc33ff", 
  "sas" = "#ffcc00"
  )

# THEMES
# Font theme
theme_nature_fonts <- function(base_size = 5) {
  theme(
    axis.text = element_text(size = base_size),
    axis.title = element_text(size = base_size),
    plot.title = element_text(size = base_size, hjust = 0.5),
    plot.subtitle = element_text(size = base_size, hjust = 0.5),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.text = element_text(size = base_size)
  )
}

# Small legend theme
theme_small_legend <- function(...) {
  theme(
    legend.key.spacing = unit(0, "cm"),
    legend.key.height = unit(0.2, "cm"),  
    legend.key.width = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
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
# Save_function
save_ggplot <- function(plot, save_path, width = 5, height = 5, with_legend = TRUE) {
  # Saves a ggplot object with or without the legend.
  # Args:
  #   plot (ggplot): The ggplot object to save.
  #   save_path (string): The file path to save the plot.
  #   width (numeric): Width of the saved plot (default: 8).
  #   height (numeric): Height of the saved plot (default: 6).
  #   with_legend (logical): Whether to keep the legend (default: TRUE).
  
  if (!with_legend) {
    plot <- plot + theme(legend.position = "none")  # Remove legend if needed
  }
  
  ggsave(filename = save_path, plot = plot, width = width, height = height)
}

# Mean variance trend plot (Voom like plot)
mean_variance_trend <- function(data, x_axis_label) {
  
  # Step 1: Compute mean and standard deviation for raw counts (no transformation)
  mean_raw <- rowMeans(data)  # Mean 
  sd_raw <- apply(data, 1, sd)  # Standard deviation 
  
  # Step 2: Create data frame for plotting
  df <- data.frame(
    mean_raw = mean_raw,
    sqrt_sd_raw = sqrt(sd_raw)  
  )

  # Step 3: Create the plot
  p <- ggplot(df, aes(x = mean_raw, y = sqrt_sd_raw)) +
    geom_point(size = 0.5) +  # Raw data points
    geom_smooth(method = "loess", color = "red", se = FALSE, linewidth = 0.5) +  # LOESS trend line
    labs(
      title = "Mean-Variance Trend",
      x = x_axis_label,  # Customizable x-axis label
      y = "Sqrt(Standard Deviation)"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme(plot.title = element_text(hjust = 0.5))  
  
  return(p)
}

# Volcano plot
volcano_plot <- function(data, logFC_thr = 1, point_size = 0.5, alpha_value = 0.5, 
                         top_features = NULL, facet_rows = NULL, facet_cols = NULL) {
  # Ensure top_features is handled correctly
  if (is.null(top_features)) {
    top_features <- character(0) 
  }
  
  # Convert character vectors to quosures for multiple facet variables
  row_facet <- if (!is.null(facet_rows)) rlang::syms(facet_rows) else NULL
  col_facet <- if (!is.null(facet_cols)) rlang::syms(facet_cols) else NULL

  p <- data |>
    ggplot(aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > logFC_thr)
      )
    ) +
    geom_point(size = point_size) +
    geom_text_repel(
      data = data |> filter(Feature %in% top_features),
      aes(label = Feature),
      size = 1.5,
      color = "black",  
      segment.color = "black",  
      min.segment.length = 0
    ) +
    facet_grid(
      rows = if (!is.null(row_facet)) vars(!!!row_facet) else NULL, 
      cols = if (!is.null(col_facet)) vars(!!!col_facet) else NULL
    ) + 
    geom_vline(
      xintercept = c(-logFC_thr, logFC_thr), 
      linetype = "dashed", 
      color = "blue",
      linewidth = (point_size/2),
      alpha = alpha_value
    ) +
    geom_hline(
      yintercept = -log10(0.05), 
      linetype = "dashed", 
      color = "blue",
      linewidth = (point_size/2),
      alpha = alpha_value
    ) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "lightgrey")
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme(legend.position = "none")

  return(p)
}

# Venn diagram
venn_diagram <- function(venn_data) {

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
    # 1. region count layer
    geom_polygon(
      data = venn_regionedge(venn_data),
      aes(X, Y, group = id),
      fill = "white"
    ) +
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
    theme(legend.position = "none")
  
  # Return the plot object
  return(p)
}

# Interaction boxplots
interactions_boxplot <- function(expression, ancestry_column, output_column, nrow, ncol) {
  
  # Generate the plot
  p <- expression %>%
    ggplot(
      aes(
        x = fct_rev(toupper(!!sym(ancestry_column))),
        y = z_score,
        fill = !!sym(output_column)
      )
    ) +
    geom_boxplot(
      width = 0.5,
      outlier.size = 0.1,
    ) +
    scale_y_continuous(
      breaks = scales::breaks_extended(n = 3)
    ) +
    facet_wrap(
      ~ Feature,   # Use `Feature` for faceting
      nrow = nrow,
      ncol = ncol,    
      scales = "free",  
      strip.position = "top"  
    ) +
    labs(
      x = "Ancestry",
      y = "z-score",
      fill = "Condition"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 60, hjust = 1)  # Rotate x-axis labels by 60 degrees
    )
  
  # Return the plot object
  return(p)
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
