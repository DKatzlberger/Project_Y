# COLORS
genetic_ancestry_colors <- c(
  "admix" = "#ff4d4d", 
  "afr"   = "#ff9900", 
  "amr"   = "#33cc33",
  "eur"   = "#3399ff", 
  "eas"   = "#cc33ff", 
  "sas"   = "#ffcc00"
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
  p <- data |>
    mutate(
      custom_x = case_when(
        !!sym(x) == levels(factor(data[[x]]))[1] ~ paste(!!sym(x), "(train_ancestry)"), 
        !!sym(x) == levels(factor(data[[x]]))[2] ~ paste(!!sym(x), "(infer_ancestry)"), 
        TRUE ~ as.character(!!sym(fill))  
      ),
      custom_fill = case_when(
        !!sym(fill) == levels(factor(data[[fill]]))[1] ~ paste(!!sym(fill), "(class_0)"), 
        !!sym(fill) == levels(factor(data[[fill]]))[2] ~ paste(!!sym(fill), "(class_1)"), 
        TRUE ~ as.character(!!sym(fill))  
      )
    ) |>
    ggplot(
      aes(
          x    = custom_x, 
          fill = custom_fill
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
      x    = paste(x, "(ancestry_column)"),
      y    = "Proportion",
      fill = paste(x, "(output_column)")
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
  p <- data |>
    mutate(
      custom_x = case_when(
        !!sym(x) == levels(factor(data[[x]]))[1] ~ paste(!!sym(x), "(train_ancestry)"), 
        !!sym(x) == levels(factor(data[[x]]))[2] ~ paste(!!sym(x), "(infer_ancestry)"), 
        TRUE ~ as.character(!!sym(fill))  
      ),
      custom_fill = case_when(
        !!sym(fill) == levels(factor(data[[fill]]))[1] ~ paste(!!sym(fill), "(class_0)"), 
        !!sym(fill) == levels(factor(data[[fill]]))[2] ~ paste(!!sym(fill), "(class_1)"), 
        TRUE ~ as.character(!!sym(fill))  
      )
    ) |>
    ggplot(
      aes(
          x    = custom_x, 
          fill = custom_fill
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
      x    = paste(x, "(ancestry_column)"),
      y    = "Count",
      fill = paste(x, "(output_column)")
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


# RESULTS
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
      x     = logFC,
      y     = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > logFC_thr)
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
    dplyr::select(Feature, patient_id, zscore) |>
    tidyr::spread(key = patient_id, value = zscore) |>
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

  
  # Ancestry colors
  genetic_ancestry_colors <- c(
    "admix" = "#ff4d4d", 
    "afr"   = "#ff9900", 
    "amr"   = "#33cc33",
    "eur"   = "#3399ff", 
    "eas"   = "#6a27a1", 
    "sas"   = "#ffcc00"
  )
  
  # Annotation
  annotations <- expression |>
    dplyr::select(patient_id, dplyr::all_of(c(output_column, ancestry_column))) |>
    dplyr::distinct() |>
    dplyr::arrange(match(patient_id, colnames(expression))) |>
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
      fill = c(genetic_ancestry_colors[train_ancestry], genetic_ancestry_colors[infer_ancestry])
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
    gp = grid::gpar(fill = genetic_ancestry_colors[train_ancestry], col = NA)
  )
  group_block_anno(
    3:4, "Ancestry", 
    gp = grid::gpar(fill = genetic_ancestry_colors[infer_ancestry], col = NA)
  )
  
  grDevices::dev.off()
  rm(heatmap, envir = .GlobalEnv)
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
