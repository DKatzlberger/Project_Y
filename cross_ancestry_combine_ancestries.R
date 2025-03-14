# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Statistics 
    library(coin)
    # Visualization
    library(patchwork)
    library(ggVennDiagram)
    }
)
# Source custom functions
source("r_utils.R")
source("figure_themes.R")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  yaml_file <- args[1] 
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_ADMIX.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
}

# Combine all ancestries for a comparison
vscratch_dir_in = file.path("data", "combined_runs")
# Construct comparison from setup file
comparison <- paste0(setup$tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))

# Save the results of the analysis
vscratch_dir_out  <- file.path("data", "combined_ancestries")
path_to_save_location <- file.path(vscratch_dir_out, comparison)
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Interactions
analysis_suffix <- "interactions"
match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")
# Extract files
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# DGE
interactions <- fload_data(match_vscratch_dir, file = "Interactions.csv")
baseline <- fload_data(match_vscratch_dir, file = "Baseline.csv")
enrichment <- fload_data(match_vscratch_dir, file = "Interactions_enrichment.csv")

# Visualize
logFC_threshold <- 1
# ---- Interactions venn (DGE) ----
gene_list <- interactions |>
  filter(adj.P.Val < 0.05, abs(logFC) > logFC_threshold) |>
  group_by(Ancestry) |>
  summarise(Features = list(Feature)) |>
  deframe()  

# Diagram
venn_data <- Venn(gene_list)
venn_data <- process_data(venn_data)
# Adjust coordinates
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

vennDiagram <- ggplot() +
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

# Save
ggsave(filename = "Interactions_venn_diagram.pdf",
       plot = vennDiagram,
       path = path_to_save_location,
       height = 2, width = 2)

# ---- Interactions enrichment ----
top_fsge <- enrichment |>
  group_by(Ancestry) |>
  arrange(desc(abs(NES))) |>
  filter(row_number() <= 15) |>
  pull(pathway)

enrichment_plot <- enrichment |>
  filter(pathway %in% top_fsge) |>
  ggplot(
    aes(
      x = Ancestry,
      y = pathway,
      color = NES,
      size = pmin(-log10(padj), 5)
    )
  ) +
  geom_point() +
  scale_size_binned(
     range = c(1, 3)    
  ) +
  scale_color_gradient2(
      high = "red", 
      mid = "white", 
      low = "blue"
  ) +
  labs(
    y = "MSigDB Hallmark 2020 gene set",
    size = "-log10(adj.P.Val)"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title.position = "top",
  )

# Save
ggsave(filename = "Interactions_fgsea.pdf",
       plot = enrichment_plot,
       path = path_to_save_location,
       height = 4, width = 4)

# ---- Baseline venn (DGE) ----
gene_list <- baseline |>
  filter(adj.P.Val < 0.05, abs(logFC) > logFC_threshold) |>
  group_by(Ancestry) |>
  summarise(Features = list(Feature)) |>
  deframe()  

# Diagram
venn_data <- Venn(gene_list)
venn_data <- process_data(venn_data)
# Adjust coordinates
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

vennDiagram <- ggplot() +
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

# Save
ggsave(filename = "Baseline_venn_diagram.pdf",
       plot = vennDiagram,
       path = path_to_save_location,
       height = 2, width = 2)



# # Analysis: Cross Ancestry -------------------------------------------------------------------------
# analysis_suffix <- "cross_ancestry"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# # Baseline -----------------------------------------------------------------------------------------
# # Summarized metric (dge)
# baseline_summarized_metric <- fload_data(
#   folders = match_vscratch_dir, 
#   file = "Baseline_summarized_metric_dge.csv"
#   )

# # Gene metric
# baseline_gene_error <- fload_data(
#   folders = match_vscratch_dir, 
#   file = "Baseline_metric_per_gene.csv"
#   )


# # Visualize ---------------------------------------------------------------------------------------------
# # Ancestry labels with sample size
# ancestry_labels <- baseline_summarized_metric |>
#   distinct(Ancestry, n_inf_ancestry) |>
#   mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")\n")) |>
#   select(Ancestry, label) |>
#   deframe()

# # Common axis
# common_y <- scale_y_continuous(
#   limits = c(0, 1.5), 
#   breaks = c(0, 0.5, 1))

# # ---- Facets prediction of average (baseline_facet_prediction_of_average_bar_plot) -----
# baseline_facet_prediction_of_average_bar_plot <- baseline_summarized_metric |>
#   mutate(Condition_with_n = paste0(Ancestry, "  ", Condition, "\n(n = ", n_condition, ")")) |>
#   ggplot(
#     aes(
#       x = fct_rev(Prediction),
#       y = mean_value
#     )
#   ) + 
#   geom_bar(
#     stat = "identity", 
#     width = 0.7
#   ) +
#   geom_errorbar(
#     aes(
#       ymin = mean_value - sd_value, 
#       ymax = mean_value + sd_value
#       ), 
#     width = 0.2, position = position_dodge(0.7)
#   ) +
#   common_y + 
#   facet_grid(
#     rows = vars(Metric),
#     col = vars(Condition_with_n)
#   ) +
#   geom_text(
#     aes(
#       x = 0.5,  # Align text to the left side of the first bar
#       y = Inf,  # Position the text at the top of the plot
#       label = ifelse(
#         is.na(p_value),
#         paste("Perm. test,", "p = NA"),  
#         paste("Perm. test,", "p = ", 
#             ifelse(p_value < 0.01, 
#                    format(p_value, digits = 3, scientific = TRUE), 
#                    format(p_value, digits = 3)))
#                    )
#       ),
#       size = 1.5,    # Adjust text size
#       vjust = 1.5,   # Align text to the top
#       hjust = 0,     # Align text to the left
#       inherit.aes = FALSE  # Don't inherit the default aesthetics
#   ) +
#   labs(
#     title = "Prediction of baseline",
#     x = "Prediction",
#     y = "Y"
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip () 

# # Save
# ggsave(filename = "Baseline_facet_prediction_of_average_bar_plot.pdf", 
#        plot = baseline_facet_prediction_of_average_bar_plot, 
#        path = path_to_save_location, 
#        width = 9, height = 3)

# # ---- Patchwork prediction of average (baseline_patchwork_prediction_of_average_bar_plot) ----
# baseline_patchwork_prediction_of_average_bar_plot_list <- list()
# for (ancestry in unique(baseline_summarized_metric$Ancestry)){
#   # Filter by ancestry
#   filtered_baseline_summarized_metric <- baseline_summarized_metric |>
#     filter(Ancestry == ancestry) |>
#     mutate(Condition_with_n = paste0(Condition, " (n = ", n_condition, ")"))


#   # Plot
#   ancestry_baseline_patchwork_prediction_of_average_bar_plot <- filtered_baseline_summarized_metric |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#       )
#     ) + 
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ), 
#       width = 0.2, position = position_dodge(0.7)
#     ) +
#     common_y + 
#     facet_grid(
#       rows = vars(Metric),
#       col = vars(Condition_with_n)
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar
#         y = Inf,  # Position the text at the top of the plot
#         label = ifelse(
#           is.na(p_value),
#           paste("p = NA"),  
#           paste("p = ", 
#               ifelse(p_value < 0.01, 
#                     format(p_value, digits = 3, scientific = TRUE), 
#                     format(p_value, digits = 3)))
#                     )
#         ),
#         size = 1.5,    # Adjust text size
#         vjust = 1.5,   # Align text to the top
#         hjust = 0,     # Align text to the left
#         inherit.aes = FALSE  # Don't inherit the default aesthetics
#     ) +
#     labs(
#       x = "Prediction",
#       y = "Y"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip () 

#   # Append to a list of plots
#   baseline_patchwork_prediction_of_average_bar_plot_list[[ancestry]] <- ancestry_baseline_patchwork_prediction_of_average_bar_plot
# }

# # Patchwork
# baseline_patchwork_prediction_of_average_bar_plot <- wrap_plots(
#   baseline_patchwork_prediction_of_average_bar_plot_list,
#   ncol = 3
#   ) 

# # Save
# ggsave(filename = "Baseline_patchwork_prediction_of_average_bar_plot.pdf", 
#        plot = baseline_patchwork_prediction_of_average_bar_plot, 
#        path = path_to_save_location, 
#        width = 6, height = 2)

# # ---- Facets gene error (baseline_facets_prediction_of_average_gene_error_boxplot) ----
# leading_genes <- baseline_gene_error |>
#   filter(Prediction == "Ancestry") |>
#   group_by(
#     Ancestry,
#     Condition,
#     Feature
#   ) |>                
#   summarize(
#     mean_RMSE = mean(RMSE, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   arrange(Ancestry, Condition, desc(mean_RMSE)) |>
#   group_by(Ancestry, Condition) |>    
#   slice_max(mean_RMSE, n = 10, with_ties = FALSE) |>
#   ungroup() |>
#   select(Ancestry, Condition, Feature)

# baseline_facets_prediction_of_average_gene_error_boxplot <- baseline_gene_error |>
#   inner_join(leading_genes, by = c("Ancestry", "Condition", "Feature")) |>
#   ggplot(
#     aes(
#       x = Feature,
#       y = RMSE
#     )
#   ) +
#   geom_boxplot(
#     outlier.size = 0.1
#   ) +
#   facet_grid(
#     cols = vars(Ancestry, fct_rev(Prediction)),
#     rows = vars(Condition),
#     scales = "free_x"
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#      axis.text.x = element_text(angle = 90, hjust = 1)
#   )

# # ---- Patchwork gene error (baseline_patchwork_prediction_of_average_gene_error_boxplot) ----
# baseline_patchwork_prediction_of_average_gene_error_boxplot_list <- list()
# for (ancestry in unique(baseline_gene_error$Ancestry)){
#   # Filter by ancestry
#   filtered_baseline_gene_error <- baseline_gene_error |>
#     filter(Ancestry == ancestry)
  
#   # Genes with highest error
#   leading_genes <- filtered_baseline_gene_error |>
#     filter(Prediction == "Ancestry") |>
#     group_by(
#       Condition,
#       Feature
#     ) |>                
#     summarize(
#       mean_RMSE = mean(RMSE, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(Condition, desc(mean_RMSE)) |>
#     group_by(Condition) |>    
#     slice_max(mean_RMSE, n = 5, with_ties = FALSE) |>
#     ungroup() |>
#     select(Condition, Feature)

#   # Plot
#   ancestry_baseline_patchwork_prediction_of_average_gene_error_boxplot <- filtered_baseline_gene_error |>
#   inner_join(leading_genes, by = c("Condition", "Feature")) |>
#   ggplot(
#     aes(
#       x = Feature,
#       y = RMSE
#     )
#   ) +
#   geom_boxplot(
#     outlier.size = 0.1
#   ) +
#   facet_grid(
#     cols = vars(fct_rev(Prediction)),
#     rows = vars(Condition),
#     scales = "free_x"
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#      axis.text.x = element_text(angle = 90, hjust = 1)
#   )

#   # Append to a list of plots
#   baseline_patchwork_prediction_of_average_gene_error_boxplot_list[[ancestry]] <- ancestry_baseline_patchwork_prediction_of_average_gene_error_boxplot
# }

# # Patchwork
# baseline_patchwork_prediction_of_average_gene_error_boxplot <- wrap_plots(
#   baseline_patchwork_prediction_of_average_gene_error_boxplot_list,
#   ncol = 3
#   ) 

# # Save
# ggsave(filename = "Baseline_patchwork_prediction_of_average_gene_error_boxplot.pdf", 
#        plot = baseline_patchwork_prediction_of_average_gene_error_boxplot, 
#        path = path_to_save_location, 
#        width = 6, height = 2)


# # Contrast ---------------------------------------------------------------------------------------------
# # Summarized metric
# contrast_summarized_metric_dge <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Summarized_metric_dge.csv"
# )
# contrast_summarized_metric_ml <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Summarized_metric_ml.csv"
# )

# # Gene metric
# contrast_gene_error <- fload(
#   folders = match_vscratch_dir,
#   file = "Contrast_metric_per_gene.csv"
# )

# # errorFC
# errorFC <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Contrast_metric_errorFC.csv"
# )

# # Visualize -------------------------------------------------------------------------------------------
# point_size = 0.5
# logFC_threshold <- 1

# # Custom labeller 
# algorithm_labels <- c(
#   "LogisticRegression" = "Logistic Regression",
#   "RandomForestClassifier" = "Random Forest")

# metric_labels <- c(
#   "ROC_AUC" = "ROC AUC",
#   "Accuracy" = "Accuracy"
# )

# # ---- Contrast patchwork correlation of logFC (contrast_patchwork_correlation_of_logFC_bar_plot) ----
# contrast_patchwork_correlation_of_logFC_bar_plot_list <- list()
# for (ancestry in unique(contrast_summarized_metric_dge$Ancestry)){
#   # Filter by ancestry
#   filtered_contrast_summarized_metric_dge <- contrast_summarized_metric_dge |>
#     filter(
#       Ancestry == ancestry,
#       Metric %in% c("Pearson", "Spearman")
#       )
  
#   # Plot
#   ancestry_contrast_patchwork_correlation_of_logFC_bar_plot <- filtered_contrast_summarized_metric_dge |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#       )
#     ) + 
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ), 
#       width = 0.2, position = position_dodge(0.7)
#     ) +
#     common_y + 
#     facet_grid(
#       cols = vars(Ancestry),
#       rows = vars(Metric),
#       labeller = labeller(Ancestry = as_labeller(ancestry_labels))
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar
#         y = Inf,  # Position the text at the top of the plot
#         label = ifelse(
#           is.na(p_value),
#           paste("p = NA"),  
#           paste("p = ", 
#               ifelse(p_value < 0.01, 
#                     format(p_value, digits = 3, scientific = TRUE), 
#                     format(p_value, digits = 3)))
#                     )
#         ),
#         size = 1.5,    # Adjust text size
#         vjust = 1.5,   # Align text to the top
#         hjust = 0,     # Align text to the left
#         inherit.aes = FALSE  # Don't inherit the default aesthetics
#     ) +
#     labs(
#       x = "Prediction",
#       y = "Y"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip () 

#   # Append to a list of plots
#   contrast_patchwork_correlation_of_logFC_bar_plot_list[[ancestry]] <- ancestry_contrast_patchwork_correlation_of_logFC_bar_plot
# }

# # Patchwork
# contrast_patchwork_correlation_of_logFC_bar_plot <-  wrap_plots(
#   contrast_patchwork_correlation_of_logFC_bar_plot_list,
#   ncol = 3
#   ) 

# # Save
# ggsave(filename = "Contrast_patchwork_correlation_of_logFC_bar_plot.pdf", 
#        plot = contrast_patchwork_correlation_of_logFC_bar_plot, 
#        path = path_to_save_location, 
#        width = 6, height = 2)

# # ---- Contrast patchwork prediction of logFC per gene (contrast_patchwork_prediction_of_logFC_bar_plot_gene_error) -----
# for (ancestry in unique(contrast_summarized_metric_dge$Ancestry)){
#   # Filter by ancestry
#   filtered_contrast_summarized_metric_dge <- contrast_summarized_metric_dge |>
#     filter(Ancestry == ancestry)

#   # Bar plot
#   ancestry_prediction_of_logFC_bar_plot <- filtered_contrast_summarized_metric_dge |>
#     filter(Metric %in% c("RMSE", "R2")) |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#         )
#     ) +
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ), 
#       width = 0.2, position = position_dodge(0.7)
#     ) +
#     common_y +
#     facet_grid(
#       rows = vars(Metric)
#     ) +
#     labs(
#       x = "Prediction",
#       y = "Y"
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar
#         y = Inf,  # Position the text at the top of the plot
#         label = ifelse(
#           is.na(p_value),
#           paste("p = NA"),  
#           paste("p = ", 
#               ifelse(p_value < 0.01, 
#                     format(p_value, digits = 3, scientific = TRUE), 
#                     format(p_value, digits = 3)))
#                     )
#         ),
#         size = 1.5,    # Adjust text size
#         vjust = 1.5,   # Align text to the top
#         hjust = 0,     # Align text to the left
#         inherit.aes = FALSE  # Don't inherit the default aesthetics
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip() 

#   # Top genes
#   genes <- contrast_gene_error |>
#     filter(Ancestry == ancestry) |>
#     filter(Prediction == "Ancestry") |>
#     group_by(
#       Comparison,
#       Feature
#     ) |>                
#     summarize(
#       mean_RMSE = mean(RMSE, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(Comparison, desc(mean_RMSE)) |>
#     group_by(Comparison) |>    
#     slice_max(mean_RMSE, n = 10) |>
#     pull(Feature) |>
#     unique()

#   # Box plot
#   ancestry_gene_error_boxplot <- contrast_gene_error |>
#     filter(Ancestry == ancestry) |>
#     filter(
#       Feature %in% genes
#     ) |>
#     ggplot(
#       aes(
#         x = Feature,
#         y = RMSE
#       )
#     ) +
#     geom_boxplot(
#       outlier.size = 0.1
#     ) +
#     facet_grid(
#       cols = vars(fct_rev(Prediction)),
#       scales = "free_x"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip() +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 10))
#     )

#   # Patchwork
#   ancestry_patchwork <- ancestry_prediction_of_logFC_bar_plot +
#     ancestry_gene_error_boxplot +
#     plot_layout(widths = c(0.8, 2)) +
#     plot_annotation(
#       title = ancestry,
#       theme = theme(plot.title = element_text(size = 5, hjust = 0.5))
#     )

#   # Save
#   plot_name <- paste0(ancestry, "_contrast_patchwork_prediction_of_logFC_bar_plot_gene_error.pdf")
#   ggsave(filename = plot_name, 
#        plot = ancestry_patchwork, 
#        path = path_to_save_location, 
#        width = 3, height = 2)
# }

# # ---- Contrast patchwork errorFC per gene (contrast_patchwork_errorFC) ----
# for (ancestry in unique(errorFC$Ancestry)){
#   # Filter by ancestry
#   filtered_errorFC <- errorFC |>
#     filter(Ancestry == ancestry)

#   # Volcano plot
#   ancestry_errorFC_volcano_plot <- filtered_errorFC |>
#     ggplot(
#       aes(
#         x = logFC,
#         y = -log10(p_adjusted),
#         color = (p_adjusted < 0.05 & abs(logFC) > logFC_threshold)
#       )
#     ) +
#     geom_point(size = point_size) +
#     geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "blue") +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
#     scale_color_manual(
#       values = c("TRUE" = "red", "FALSE" = "grey"),
#     ) +
#     labs(
#       y = "-log10(adj.P.Val)"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip() +
#     theme(legend.position = "none") 
  
#   # Bar plot
#   ancestry_errorFC_bar_plot <- filtered_errorFC |>
#     filter(p_adjusted < 0.05 & abs(logFC) > logFC_threshold) |>
#     slice_max(logFC, n = 10) |>
#     ggplot(
#       aes(
#         x = Feature,
#         y = logFC
#       )
#     ) +
#     geom_col(
#       width = 0.7
#     ) +
#     labs(
#       x = "Feature",
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 10))
#     )

#     # Patchwork
#   ancestry_patchwork <- ancestry_errorFC_volcano_plot +
#     ancestry_errorFC_bar_plot +
#     plot_layout(widths = c(0.8, 2)) +
#     plot_annotation(
#       title = ancestry,
#       theme = theme(plot.title = element_text(size = 5, hjust = 0.5))
#     )

#   # Save
#   plot_name <- paste0(ancestry, "_contrast_patchwork_errorFC.pdf")
#   ggsave(filename = plot_name, 
#        plot = ancestry_patchwork, 
#        path = path_to_save_location, 
#        width = 3, height = 2)
# }

# # ---- Contrast patchwork prediction of phenotype (contrast_patchwork_prediction_of_phenotype_bar_plot) -----
# contrast_patchwork_prediction_of_phenotype_bar_plot_list <- list()
# for (ancestry in unique(contrast_summarized_metric_ml$Ancestry)){
  
#   # Filter by ancestry
#   filtered_contrast_summarized_metric_ml <- contrast_summarized_metric_ml |>
#     filter(
#       Ancestry == ancestry
#       )
  
#   # Plot
#   ancestry_contrast_patchwork_prediction_of_phenotype_bar_plot <- filtered_contrast_summarized_metric_ml |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#       )
#     ) + 
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ), 
#       width = 0.2, position = position_dodge(0.7)
#     ) +
#     common_y + 
#     facet_grid(
#       cols = vars(Algorithm),
#       rows = vars(Metric),
#       labeller = labeller(
#         Algorithm = as_labeller(algorithm_labels),
#         Metric = as_labeller(metric_labels)
#         )
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar
#         y = Inf,  # Position the text at the top of the plot
#         label = ifelse(
#           is.na(p_value),
#           paste("p = NA"),  
#           paste("p = ", 
#               ifelse(p_value < 0.01, 
#                     format(p_value, digits = 3, scientific = TRUE), 
#                     format(p_value, digits = 3)))
#                     )
#         ),
#         size = 1.5,    # Adjust text size
#         vjust = 1.5,   # Align text to the top
#         hjust = 0,     # Align text to the left
#         inherit.aes = FALSE  # Don't inherit the default aesthetics
#     ) +
#     labs(
#       x = "Prediction",
#       y = "Y"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip ()

#   # Append to a list of plots
#   contrast_patchwork_prediction_of_phenotype_bar_plot_list[[ancestry]] <- ancestry_contrast_patchwork_prediction_of_phenotype_bar_plot
# }

# # Patchwork
# contrast_patchwork_prediction_of_phenotype_bar_plot <-  wrap_plots(
#   contrast_patchwork_prediction_of_phenotype_bar_plot_list,
#   ncol = 3
#   ) 

# # Save
# ggsave(filename = "Contrast_patchwork_prediction_of_phenotype_bar_plot.pdf", 
#        plot = contrast_patchwork_prediction_of_phenotype_bar_plot, 
#        path = path_to_save_location, 
#        width = 6, height = 2)



# # Analysis: Interactions -------------------------------------------------------------------------------
# analysis_suffix <- "interactions"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)


# # Baseline/interactions -------------------------------------------------------------------------------------------
# baseline <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Baseline.csv"
# )

# interactions <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Interactions.csv"
# )

# # Visualize ---------------------------------------------------------------------------------------
# point_size = 0.5
# logFC_threshold <- 1

# # ---- Baseline facet volcano plot (baseline_facet_volcano_plot) -----
# baseline_facet_volcano_plot <- baseline |>
#   ggplot(
#     aes(
#       x = logFC,
#       y = -log10(adj.P.Val),
#       color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
#     )
#   ) +
#   geom_point(size = point_size) +
#   facet_grid(
#     cols = vars(Ancestry),
#     rows = vars(Condition),
#     labeller = labeller(Ancestry = as_labeller(ancestry_labels))
#   ) + 
#   geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
#   scale_color_manual(
#     values = c("TRUE" = "red", "FALSE" = "grey"),
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#     legend.position = "none"
#   )

# # Save
# ggsave(filename = "Baseline_facet_volcano_plot.pdf", 
#        plot = baseline_facet_volcano_plot, 
#        path = path_to_save_location, 
#        width = 6, height = 2)

# # ---- Interactions facet volcano plot (interactions_facet_volcano_plot) ----
# interactions_facet_volcano_plot <- interactions |>
#   ggplot(
#     aes(
#       x = logFC,
#       y = -log10(adj.P.Val),
#       color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
#     )
#   ) +
#   geom_point(size = point_size) +
#   facet_grid(
#     cols = vars(Ancestry),
#     labeller = labeller(Ancestry = as_labeller(ancestry_labels))
#   ) + 
#   geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
#   scale_color_manual(
#     values = c("TRUE" = "red", "FALSE" = "grey"),
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#     legend.position = "none"
#   )

# # Save
# ggsave(filename = "Interactions_facet_volcano_plot.pdf", 
#        plot = interactions_facet_volcano_plot, 
#        path = path_to_save_location, 
#        width = 3, height = 2)

# # ---- Interactions patchwork volcano plot + gene logFC (interactions_patchwork_volcano_plot_gene_error) ----
# for (ancestry in unique(interactions$Ancestry)){
#   # Filter by ancestry
#   filtered_interactions <- interactions |>
#     filter(Ancestry == ancestry)
  
#   # Volcano plot
#   ancestry_volcano_plot <- filtered_interactions |>
#     ggplot(
#       aes(
#         x = logFC,
#         y = -log10(adj.P.Val),
#         color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
#       )
#     ) +
#     geom_point(size = point_size) +
#     geom_vline(
#       xintercept = c(-logFC_threshold, logFC_threshold), 
#       linetype = "dashed", 
#       color = "blue",
#       # linewidth = point_size
#     ) +
#     geom_hline(
#       yintercept = -log10(0.05), 
#       linetype = "dashed", 
#       color = "blue",
#       # linewidth = point_size
#     ) +
#     scale_color_manual(
#       values = c("TRUE" = "red", "FALSE" = "grey"),
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme_white_strip() +
#     theme(
#       legend.position = "none"
#     )
  
#   # Significant interactions
#   ancestry_significant_logFC_plot <- filtered_interactions |>
#     filter(adj.P.Val < 0.05 & abs(logFC) > logFC_threshold) |>
#     ggplot(
#       aes(
#         x = Feature,
#         y = logFC
#       )
#     ) +
#     geom_col(width = 0.7) +
#     labs(
#       x = "Feature",
#       y = "logFC"
#     ) +
#     theme_nature_fonts() +
#     theme_white_background() +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 10))
#     )

#   # Patchwork
#   ancestry_patchwork <- ancestry_volcano_plot + 
#     ancestry_significant_logFC_plot +
#     plot_layout(widths = c(0.8, 2)) +
#     plot_annotation(
#       title = ancestry,
#       theme = theme(plot.title = element_text(size = 5, hjust = 0.5))
#     )
  
#   # Save
#   plot_name <- paste0(ancestry, "_interactions_patchwork_volcano_plot_gene_error.pdf")
#   ggsave(filename = plot_name, 
#        plot = ancestry_patchwork, 
#        path = path_to_save_location, 
#        width = 3, height = 2)
# }


# # Analysis: Subsetting -------------------------------------------------------------------------------
# analysis_suffix <- "subsetting"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# # Load data
# subsetting <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Combined_metric.csv"
# )

# # Visualize ---------------------------------------------------------------------------------------------
# # ---- Performance by sample size (performance_by_sample_line_plot) ----
# subsetting <- subsetting |>
#   mutate(
#     Metric_type = recode(Metric_type,
#                          "ROC_AUC (LogisticRegression)" = "ROC AUC (Logistic Regression)",
#                          "ROC_AUC (RandomForestClassifier)" = "ROC AUC (Random Forest)",
#                          "Spearman (limma)" = "Spearman (Limma)",
#                          "Pearson (limma)" = "Pearson (Limma)"),
#     Metric_type = factor(Metric_type, levels = c(
#       "ROC AUC (Logistic Regression)", 
#       "ROC AUC (Random Forest)", 
#       "Spearman (Limma)",
#       "Pearson (Limma)"
#     ))
#   )

# performance_by_sample_line_plot <- subsetting |>
#   ggplot(
#     aes(
#       x = n_test_ancestry,
#       y = mean_value,
#       color = Metric_type
#     )
#   ) +
#   geom_point(size = point_size) +
#   geom_line(linewidth = point_size) +
#   geom_errorbar(
#     aes(
#       ymin = mean_value - sd_value, 
#       ymax = mean_value + sd_value
#       ), 
#     width = 0.1
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1.1),
#     breaks = c(0.0, 0.5, 1.0)
#   ) +
#   labs(
#     x = "n (Testset)",
#     y = "Y",
#     color = "Method"
#   ) +
#   theme_nature_fonts() +
#   theme_small_legend() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = c(1, 0.1),  
#     legend.direction = "vertical",
#     legend.justification = c(1, 0)
#   ) 

# # Save
# ggsave(filename = "Performance_by_sample_line_plot.pdf", 
#        plot = performance_by_sample_line_plot, 
#        path = path_to_save_location, 
#        width = 3, height = 2)


# # Analysis: Robustness -------------------------------------------------------------------------------
# analysis_suffix <- "robustness"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# # Load data
# robustness_summarized_metric_dge <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Summarized_metric_dge.csv"
# )

# robustness_summarized_metric_ml <- fload_data(
#   folders = match_vscratch_dir,
#   file = "Summarized_metric_ml.csv"
# )

# # Visualize -----------------------------------------------------------------------------------------------
# algorithm_labels <- c(
#   "LogisticRegression" = "LR",
#   "RandomForestClassifier" = "RF")

# # ---- Robustness patchwork correaltion of logFC (robustness_patchwork_correlation_of_logFC_bar_plot) -----
# for (ancestry in unique(robustness_summarized_metric_dge$Ancestry)){
#   prop_robustness_patchwork_correlation_of_logFC_bar_plot_list <- list()
#   for (prop in unique(robustness_summarized_metric_dge$Proportion)){
#     # Filter for ancestry and proportion
#     filtered_robustness_summarized_metric_dge <- robustness_summarized_metric_dge |>
#       filter(
#         Ancestry == ancestry,
#         Proportion == prop)
    
#     # Label for sample size of proportion
#     n_inf_label <- filtered_robustness_summarized_metric_dge |>
#       ungroup() |>  
#       distinct(Ancestry, n_inf_ancestry, Proportion) |>
#       mutate(label = paste0("n = ", n_inf_ancestry, " (", Proportion, ")")) |>
#       select(Ancestry, label) |>
#       deframe()

#     # Plot
#     ancestry_prop_robustness_patchwork_correlation_of_logFC_bar_plot <- filtered_robustness_summarized_metric_dge |>
#       ggplot(
#         aes(
#           x = fct_rev(Prediction),
#           y = mean_value
#         )
#       ) + 
#       geom_bar(
#         stat = "identity", 
#         width = 0.7
#       ) +
#       geom_errorbar(
#         aes(
#           ymin = mean_value - sd_value, 
#           ymax = mean_value + sd_value
#           ), 
#         width = 0.2, position = position_dodge(0.7)
#       ) +
#       common_y + 
#       facet_grid(
#         rows = vars(Metric),
#         col = vars(Ancestry),
#         labeller = labeller(Ancestry = as_labeller(n_inf_label))
#       ) +
#       geom_text(
#         aes(
#           x = 0.5,  # Align text to the left side of the first bar
#           y = Inf,  # Position the text at the top of the plot
#           label = ifelse(
#             is.na(p_value),
#             paste("p = NA"),  
#             paste("p = ", 
#                 ifelse(p_value < 0.01, 
#                       format(p_value, digits = 3, scientific = TRUE), 
#                       format(p_value, digits = 3)))
#                       )
#           ),
#           size = 1.5,    # Adjust text size
#           vjust = 1.5,   # Align text to the top
#           hjust = 0,     # Align text to the left
#           inherit.aes = FALSE  # Don't inherit the default aesthetics
#       ) +
#       labs(
#         x = "Prediction",
#         y = "Y"
#       ) +
#       theme_nature_fonts() +
#       theme_white_background() +
#       theme_white_strip () +
#       theme(plot.margin = margin(0, 0, 0, 0))

#       # Append to a list of plots
#       prop_robustness_patchwork_correlation_of_logFC_bar_plot_list[[paste0("prop_", prop)]] <- ancestry_prop_robustness_patchwork_correlation_of_logFC_bar_plot
#   }

#   # 1st Patchwork layer
#   ancestry_robustness_patchwork_correlation_of_logFC_bar_plot <- wrap_plots(
#     prop_robustness_patchwork_correlation_of_logFC_bar_plot_list, 
#     ncol = 2
#     ) +
#     plot_annotation(
#       title = ancestry,
#       theme = theme(plot.title = element_text(size = 5, hjust = 0.5))
#     )
  
#   # Save 
#   plot_name <- paste0(ancestry, "_robustness_patchwork_correlation_of_logFC_bar_plot.pdf")
#   ggsave(filename = plot_name, 
#        plot = ancestry_robustness_patchwork_correlation_of_logFC_bar_plot, 
#        path = path_to_save_location, 
#        width = 3, height = 3)
# }

# # ---- Robustness patchwork prediction of phenotype (robustness_patchwork_prediction_of_phenotype_bar_plot) -----
# for (ancestry in unique(robustness_summarized_metric_ml$Ancestry)){
#   prop_robustness_patchwork_prediction_of_phenotype_bar_plot_list <- list()
#   for (prop in unique(robustness_summarized_metric_ml$Proportion)){
#     # Filter for ancestry and proportion
#     filtered_robustness_summarized_metric_ml <- robustness_summarized_metric_ml |>
#       filter(
#         Ancestry == ancestry,
#         Proportion == prop)
    
#     # Label for sample size of proportion
#     n_inf_label <- filtered_robustness_summarized_metric_ml |>
#       ungroup() |>  
#       distinct(Ancestry, n_inf_ancestry, Proportion) |>
#       mutate(label = paste0("n = ", n_inf_ancestry, " (", Proportion, ")")) |>
#       select(Ancestry, label) |>
#       deframe()

#     # Plot
#     ancestry_prop_robustness_patchwork_prediction_of_phenotype_bar_plot <- filtered_robustness_summarized_metric_ml |>
#       ggplot(
#         aes(
#           x = fct_rev(Prediction),
#           y = mean_value
#         )
#       ) + 
#       geom_bar(
#         stat = "identity", 
#         width = 0.7
#       ) +
#       geom_errorbar(
#         aes(
#           ymin = mean_value - sd_value, 
#           ymax = mean_value + sd_value
#           ), 
#         width = 0.2, position = position_dodge(0.7)
#       ) +
#       common_y + 
#       facet_grid(
#         cols = vars(Ancestry),
#         rows = vars(Algorithm),
#         labeller = labeller(
#           Ancestry = as_labeller(n_inf_label),
#           Algorithm = as_labeller(algorithm_labels)
#           )
#       ) +
#       geom_text(
#         aes(
#           x = 0.5,  # Align text to the left side of the first bar
#           y = Inf,  # Position the text at the top of the plot
#           label = ifelse(
#             is.na(p_value),
#             paste("p = NA"),  
#             paste("p = ", 
#                 ifelse(p_value < 0.01, 
#                       format(p_value, digits = 3, scientific = TRUE), 
#                       format(p_value, digits = 3)))
#                       )
#           ),
#           size = 1.5,    # Adjust text size
#           vjust = 1.5,   # Align text to the top
#           hjust = 0,     # Align text to the left
#           inherit.aes = FALSE  # Don't inherit the default aesthetics
#       ) +
#       labs(
#         x = "Prediction",
#         y = "Y"
#       ) +
#       theme_nature_fonts() +
#       theme_white_background() +
#       theme_white_strip () +
#       theme(plot.margin = margin(0, 0, 0, 0))

#       # Append to a list of plots
#       prop_robustness_patchwork_prediction_of_phenotype_bar_plot_list[[paste0("prop_", prop)]] <- ancestry_prop_robustness_patchwork_prediction_of_phenotype_bar_plot
#   }
#   # 1st Patchwork layer
#   ancestry_robustness_patchwork_prediction_of_phenotype_bar_plot <- wrap_plots(
#     prop_robustness_patchwork_prediction_of_phenotype_bar_plot_list, 
#     ncol = 2) +
#     plot_annotation(
#       title = ancestry,
#       theme = theme(plot.title = element_text(size = 5, hjust = 0.5))
#     )
  
#   # Save 
#   plot_name <- paste0(ancestry, "_robustness_patchwork_prediction_of_phenotype_bar_plot.pdf")
#   ggsave(filename = plot_name, 
#        plot = ancestry_robustness_patchwork_prediction_of_phenotype_bar_plot, 
#        path = path_to_save_location, 
#        width = 3, height = 3)
# }


# # Functional analysis ----------------------------------------------------------------------------------------
# # Load enriched limma interactions
# analysis_suffix <- "interactions"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# # Load data
# interactions_enrichment <- fload_data(
#   folders = match_vscratch_dir, 
#   file = "Interactions_enrichment.csv"
#   )
# interactions_enrichment[, ranked_by := "Interaction logFC"]

# # Load enriched errorFC
# analysis_suffix <- "cross_ancestry"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# # Load data
# contrast_enrichment <- fload_data(
#   folders = match_vscratch_dir, 
#   file = "Contrast_enrichment.csv"
#   ) 
# contrast_enrichment[, ranked_by := "Error logFC"]

# # Combine
# combined_enrichment <- bind_rows(interactions_enrichment, contrast_enrichment)

# # Visualize -------------------------------------------------------------------------------------
# # ---- Dotplot of enrichments (factes_fgsea_enrichment_dotplot) -----
# top_pathways <- interactions_enrichment |>
#   arrange(desc(abs(NES))) |>
#   filter(row_number() <= 30) |>
#   pull(pathway)
  
# # Dotplot
# factes_fgsea_enrichment_dotplot <- combined_enrichment |>
#   filter(pathway %in% top_pathways) |>
#   ggplot(
#     aes(
#       x = ranked_by,
#       y = pathway,
#       color = NES,
#       size = pmin(-log10(padj), 5)
#     )
#   ) +
#   geom_point() +
#   facet_grid(
#     cols = vars(Ancestry)
#   ) +
#   scale_size_binned(
#     range = c(1, 3)    
#   ) +
#   scale_color_gradient2(
#       high = "red", 
#       mid = "white", 
#       low = "blue"
#   ) +
#   labs(
#     x = "Gene ranking",
#     y = "MSigDB Hallmark 2020 gene set",
#     size = "-log10(adj.P.Val)"
#   ) +
#   theme_nature_fonts() +
#   theme_small_legend() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.title.position = "top"
#   )

# # Save
# ggsave(filename = "Facets_fgsea_enrichment_dotplot.pdf", 
#        plot = factes_fgsea_enrichment_dotplot, 
#        path = path_to_save_location, 
#        width = 6, height = 4)




# all_combined_plots <- list()
# for (ancestry in ancestries){
#   ggplot_list_dge = list()
#   for (prop in proportions) {
#     # Filter by proportion
#     filtered_dge <- summarized_metric_dge |> filter(
#       Proportion == prop,
#       Ancestry == ancestry
#       )

#     # Create lable to showcase number of samples
#     n_inf_label <- filtered_dge |>
#       ungroup() |>  
#       distinct(Ancestry, n_inf_ancestry) |>
#       mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#       select(Ancestry, label) |>
#       deframe()
    
#     # Correlation axis from 0 to 1
#     common_y <- scale_y_continuous(
#       limits = c(0, 1.2), 
#       breaks = c(0, 0.5, 1))
    
#     # Create plot
#     prop_plot <- filtered_dge |>
#       ggplot(
#         aes(
#           x = fct_rev(Prediction),
#           y = mean_value
#           )
#       ) +
#       geom_bar(
#         stat = "identity", 
#         width = 0.7
#       ) +
#       geom_errorbar(
#         aes(
#           ymin = mean_value - sd_value, 
#           ymax = mean_value + sd_value
#           ),
#         width = 0.2, 
#         position = position_dodge(0.7)
#       ) +
#       common_y +
#       facet_grid(
#         rows = vars(Metric),
#         col = vars(Ancestry),
#         labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                             Metric = label_value)
#       ) +
#       labs(
#           # title = paste(gsub("_", " ", comparison)),
#           x = "Prediction",
#           y = "Correlation coefficient"
#       ) +
#       geom_text(
#         aes(
#           x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
#           y = Inf,  # Position the text at the top of the plot
#           label = paste("Perm. test,", "p.adj =", format(p_adjusted, digits = 3))
#           ),
#         size = 3,    
#         vjust = 1.5,   # Align text to the top (since we're using Inf for y)
#         hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
#         inherit.aes = FALSE  
#     ) 

#     # Add the plot to the list with a unique name
#     ggplot_list_dge[[paste0("prop_", prop)]] <- prop_plot
#   }
#   # Combine with patchwork 
#   combined_dge_bar_plot <- wrap_plots(ggplot_list_dge, ncol = 2)
#   # Combine ancestry wise
#   all_combined_plots[[ancestry]] <- combined_dge_bar_plot
# }




# # Combining metric csv files 
# # Each folder is one seed
# summarized_metric_dge <- data.frame()
# summarized_metric_ml <- data.frame()
# for (folder in match_vscratch_dir){
#     dge_file <- file.path(folder, "Summarized_metric_dge.csv")
#     ml_file <- file.path(folder, "Summarized_metric_ml.csv")

#     # Check if files exists
#     if (file.exists(dge_file)) {
#       # Read data
#       dge_data <- fread(dge_file)
#       # Append data
#       summarized_metric_dge <- bind_rows(summarized_metric_dge, dge_data)
#     } else {
#       warning("File does not exist: ", dge_file)
#     }

#     # Check if file exists
#     if (file.exists(ml_file)) {
#       # Read data
#       ml_data <- fread(ml_file)
#       # Append data
#       summarized_metric_ml <- bind_rows(summarized_metric_ml, ml_data)
#     } else {
#       warning("File does not exist: ", ml_file)
#     }
# }

# # Pvalue correction
# # dge
# unique_p_values_dge <- summarized_metric_dge |> pull(p_value) |> unique()
# # Benjamini-Hochberg correction
# p_adjusted_bh_dge <- p.adjust(unique_p_values_dge, method = "BH")
# # Create a named vector to map the adjusted p-values back to the data frame
# p_adjusted_named_dge <- setNames(p_adjusted_bh_dge, unique_p_values_dge)
# # Add 
# summarized_metric_dge <- summarized_metric_dge |>
#    mutate(p_adjusted = p_adjusted_named_dge[as.character(p_value)])

# # ml
# unique_p_values_ml <- summarized_metric_ml |> pull(p_value) |> unique()
# # Benjamini-Hochberg correction
# p_adjusted_bh_ml <- p.adjust(unique_p_values_ml, method = "BH")
# # Create a named vector to map the adjusted p-values back to the data frame
# p_adjusted_named_ml <- setNames(p_adjusted_bh_ml, unique_p_values_ml)
# # Add 
# summarized_metric_ml <- summarized_metric_ml |>
#    mutate(p_adjusted = p_adjusted_named_ml[as.character(p_value)])

# # Save the metric data frames
# fwrite(summarized_metric_dge, file.path(path_to_save_location, "Summarized_metric_dge.csv"))
# fwrite(summarized_metric_ml, file.path(path_to_save_location, "Summarized_metric_ml.csv"))


# # Visualization:
# # Bar plots
# ancestries <- unique(summarized_metric_dge$Ancestry)
# # dge
# ggplot_list_dge = list()
# for (ancestry in ancestries){
#   filtered_metric <- filter(summarized_metric_dge, Ancestry == ancestry)
#   # Create lable to showcase number of samples
#   n_inf_label <- filtered_metric |>
#     ungroup() |>  
#     distinct(Ancestry, n_inf_ancestry) |>
#     mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#     select(Ancestry, label) |>
#     deframe()

#   # Correlation axis from 0 to 1
#   common_y <- scale_y_continuous(
#     limits = c(0, 1.2), 
#     breaks = c(0, 0.5, 1))

#   plot <- filtered_metric |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#         )
#     ) +
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ),
#       width = 0.2, 
#       position = position_dodge(0.7)
#     ) +
#     common_y +
#     facet_grid(
#       rows = vars(Metric),
#       col = vars(Ancestry),
#       labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                           Metric = label_value)
#     ) +
#     labs(
#         # title = paste(gsub("_", " ", comparison)),
#         x = "Prediction",
#         y = "Correlation coefficient"
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
#         y = Inf,  # Position the text at the top of the plot
#         label = paste("Perm. test,", "p.adj =", format(p_adjusted, digits = 3))
#         ),
#       size = 3,    
#       vjust = 1.5,   # Align text to the top (since we're using Inf for y)
#       hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
#       inherit.aes = FALSE  
#   ) 

#   # Add the plot to the list with a unique name
#   ggplot_list_dge[[paste0("ancestry_", ancestry)]] <- plot
# }

# # Combine with patchwork 
# combined_dge_bar_plot <- wrap_plots(ggplot_list_dge, ncol = 3) +
#   plot_annotation(
#     title = "Generalizability of logFC across ancestries"
#   ) &
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )

# # Save the image
# ggsave(filename = "Plot_bar_dge.pdf", 
#        plot = combined_dge_bar_plot, 
#        path = path_to_save_location, 
#        width = 8, height = 4)

# # Density plots
# ancestries <- unique(metric_dge$Ancestry)
# # dge
# ggplot_list_dge = list()
# for (ancestry in ancestries){
#   filtered_metric <- filter(metric_dge, Ancestry == ancestry)
#   # Create lable to showcase number of samples
#   n_inf_label <- filtered_metric |>
#     ungroup() |>  
#     distinct(Ancestry, n_inf_ancestry) |>
#     mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#     select(Ancestry, label) |>
#     deframe()
  
#   # Correlation axis from 0 to 1
#   common_x <- scale_x_continuous(
#     limits = c(0.5, 1.2), 
#     breaks = c(0.5, 1))

#   # Plot
#   plot <- filtered_metric |>
#     pivot_longer(
#     cols = c(Pearson, Spearman), 
#     names_to = "Metric", 
#     values_to = "Value"
#     ) |>
#     ggplot(
#       aes(
#         x = Value,
#         fill = fct_rev(Prediction)
#       )
#     ) +
#     geom_density(
#       alpha = 0.5
#     ) +
#     common_x +
#     facet_grid(
#       rows = vars(Metric),
#       col = vars(Ancestry),
#       labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                           Metric = label_value),
#       scales = "free"
#     ) +
#     labs(
#       x = "Correlation coefficient",
#       y = "Density",
#       fill = "Prediction"
#     ) +
#     theme(
#       legend.position = "bottom",
#       legend.direction = "horizontal"
#     )
  
#   # Add the plot to the list with a unique name
#   ggplot_list_dge[[paste0("ancestry_", ancestry)]] <- plot
# }
# # Combine with patchwork 
# combined_dge_density_plot <- wrap_plots(ggplot_list_dge, ncol = 3) + 
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")

# # Save the plot
# ggsave(filename = "Plot_density_dge.pdf", 
#        plot = combined_dge_density_plot, 
#        path = path_to_save_location, 
#        width = 8, height = 4)

# # Density and Bar
# dge_bar_density <- combined_dge_bar_plot / combined_dge_density_plot

# # Save the plot
# ggsave(filename = "Plot_bar_density_dge.pdf", 
#        plot = dge_bar_density, 
#        path = path_to_save_location, 
#        width = 7.2, height = 6)



# # ml 
# ggplot_list_ml = list()
# for (ancestry in ancestries) {
#   filtered_metric <- filter(summarized_metric_ml, Ancestry == ancestry)

#   # Create lable to showcase number of samples
#   n_inf_label <- filtered_metric |>
#     ungroup() |>  
#     distinct(Ancestry, n_inf_ancestry) |>
#     mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#     select(Ancestry, label) |>
#     deframe()
  
#   # Correlation axis from 0 to 1
#   common_y <- scale_y_continuous(
#     limits = c(0, 1.2), 
#     breaks = c(0, 0.5, 1))
  
#   # Create plot
#   plot <- filtered_metric |>
#     ggplot(
#       aes(
#         x = fct_rev(Prediction),
#         y = mean_value
#         )
#     ) +
#     geom_bar(
#       stat = "identity", 
#       width = 0.7
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = mean_value - sd_value, 
#         ymax = mean_value + sd_value
#         ),
#       width = 0.2, 
#       position = position_dodge(0.7)
#     ) +
#     common_y +
#     facet_grid(
#       rows = vars(Algorithm),
#       col = vars(Ancestry),
#       labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                           Algorithm = label_value)
#     ) +
#     labs(
#         # title = paste(gsub("_", " ", comparison)),
#         x = "Prediction",
#         y = "ROC AUC"
#     ) +
#     geom_text(
#       aes(
#         x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
#         y = Inf,  # Position the text at the top of the plot
#         label = paste("Perm. test,", "p.adj =", format(p_adjusted, digits = 3))
#         ),
#       size = 3,    
#       vjust = 1.5,   # Align text to the top (since we're using Inf for y)
#       hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
#       inherit.aes = FALSE  
#   ) 

#   # Add the plot to the list with a unique name
#   ggplot_list_ml[[paste0("ancestry_", ancestry)]] <- plot
# }

# # Combine with patchwork 
# combined_ml_bar_plot <- wrap_plots(ggplot_list_ml, ncol = 3)
# # Save the image
# ggsave(filename = "Plot_bar_ml.pdf", 
#        plot = combined_ml_bar_plot, 
#        path = path_to_save_location, 
#        width = 8, height = 4)


# # Density
# # ml
# ggplot_list_ml = list()
# for (ancestry in ancestries){
#   filtered_metric <- filter(metric_ml, Ancestry == ancestry)
#   # Create lable to showcase number of samples
#   n_inf_label <- filtered_metric |>
#     ungroup() |>  
#     distinct(Ancestry, n_inf_ancestry) |>
#     mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#     select(Ancestry, label) |>
#     deframe()
  
#   # Correlation axis from 0 to 1
#   common_x <- scale_x_continuous(
#     limits = c(0.5, 1.2), 
#     breaks = c(0.5, 1))

#   # Plot
#   plot <- filtered_metric |>
#     pivot_longer(
#       cols = c(LogisticRegression, RandomForestClassifier), 
#       names_to = "Algorithm", 
#       values_to = "Value"
#     ) |>
#     ggplot(
#       aes(
#         x = Value,
#         fill = fct_rev(Prediction)
#       )
#     ) +
#     geom_density(
#       alpha = 0.5
#     ) +
#     common_x +
#     facet_grid(
#       rows = vars(Algorithm),
#       col = vars(Ancestry),
#       labeller = labeller(Ancestry = as_labeller(n_inf_label)),
#       scales = "free"
#     ) +
#     labs(
#       x = "ROC AUC",
#       y = "Density",
#       fill = "Prediction"
#     ) +
#     theme(
#       legend.position = "bottom",
#       legend.direction = "horizontal"
#     )
  
#   # Add the plot to the list with a unique name
#   ggplot_list_ml[[paste0("ancestry_", ancestry)]] <- plot
# }
# # Combine with patchwork 
# combined_ml_density_plot <- wrap_plots(ggplot_list_ml, ncol = 3) + 
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")

# # Save the plot
# ggsave(filename = "Plot_density_ml.pdf", 
#        plot = combined_ml_density_plot, 
#        path = path_to_save_location, 
#        width = 8, height = 4)


# # Density and Bar
# ml_bar_density <- combined_ml_bar_plot / combined_ml_density_plot

# # Save the plot
# ggsave(filename = "Plot_bar_density_ml.pdf", 
#        plot = ml_bar_density, 
#        path = path_to_save_location, 
#        width = 7.2, height = 6)


# # Combine bar plots
# dge_bar_ml_bar <- combined_dge_bar_plot / combined_ml_bar_plot

# # Save the plot
# ggsave(filename = "Plot_bar_dge_bar_ml.pdf", 
#        plot = dge_bar_ml_bar, 
#        path = path_to_save_location, 
#        width = 7.3, height = 7.3)






# # --------------------------------------------------------------------------------------------
# # Analysis: Subsetting
# analysis_suffix <- "subsetting"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# subsetting <- data.frame()
# for (folder in match_vscratch_dir){
#     subsetting_file <- file.path(folder, "Combined_metric.csv")
  
#     # Check if files exists
#     if (file.exists(subsetting_file)) {
#       # Read data
#       subsetting_data <- fread(subsetting_file)
#       # Append data
#       subsetting <- bind_rows(subsetting, subsetting_data)
#     } else {
#       warning("File does not exist: ", interaction_file)
#     }
# }

# subsetting <- subsetting |>
#   mutate(
#     Metric_type = factor(Metric_type, 
#                          levels = c(
#                            "ROC_AUC (LogisticRegression)", 
#                            "ROC_AUC (RandomForestClassifier)", 
#                            "Spearman (limma)",
#                            "Pearson (limma)"
#                          )
#                         )
#   )

# performance_plot <- subsetting |>
#   ggplot(
#     aes(
#       x = n_test_ancestry,
#       y = mean_value,
#       color = Metric_type
#     )
#   ) +
#   geom_point(size = 2) +
#   geom_line() +
#   geom_errorbar(
#     aes(
#       ymin = mean_value - sd_value, 
#       ymax = mean_value + sd_value
#       ), 
#     width = 0.1
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1.1),
#     breaks = c(0.0, 0.5, 1.0)
#   ) +
#   labs(
#     title = "Performance by sample size",
#     x = "Test sample size (EUR)",
#     y = "Mean value",
#     color = "Method"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )



# # --------------------------------------------------------------------------------------------
# # Analysis: Robustness
# analysis_suffix <- "robustness"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)


# # Load robustness files
# summarized_metric_dge <- data.frame()
# summarized_metric_ml <- data.frame()
# for (folder in match_vscratch_dir){
#     dge_file <- file.path(folder, "Summarized_metric_dge.csv")
#     ml_file <- file.path(folder, "Summarized_metric_ml.csv")

#     # Check if files exists
#     if (file.exists(dge_file)) {
#       # Read data
#       dge_data <- fread(dge_file)
#       # Append data
#       summarized_metric_dge <- bind_rows(summarized_metric_dge, dge_data)
#     } else {
#       warning("File does not exist: ", dge_file)
#     }

#     # Check if file exists
#     if (file.exists(ml_file)) {
#       # Read data
#       ml_data <- fread(ml_file)
#       # Append data
#       summarized_metric_ml <- bind_rows(summarized_metric_ml, ml_data)
#     } else {
#       warning("File does not exist: ", ml_file)
#     }
# }

# proportions <- summarized_metric_dge |> pull(Proportion) |> unique()

# all_combined_plots <- list()
# for (ancestry in ancestries){
#   ggplot_list_dge = list()
#   for (prop in proportions) {
#     # Filter by proportion
#     filtered_dge <- summarized_metric_dge |> filter(
#       Proportion == prop,
#       Ancestry == ancestry
#       )

#     # Create lable to showcase number of samples
#     n_inf_label <- filtered_dge |>
#       ungroup() |>  
#       distinct(Ancestry, n_inf_ancestry) |>
#       mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#       select(Ancestry, label) |>
#       deframe()
    
#     # Correlation axis from 0 to 1
#     common_y <- scale_y_continuous(
#       limits = c(0, 1.2), 
#       breaks = c(0, 0.5, 1))
    
#     # Create plot
#     prop_plot <- filtered_dge |>
#       ggplot(
#         aes(
#           x = fct_rev(Prediction),
#           y = mean_value
#           )
#       ) +
#       geom_bar(
#         stat = "identity", 
#         width = 0.7
#       ) +
#       geom_errorbar(
#         aes(
#           ymin = mean_value - sd_value, 
#           ymax = mean_value + sd_value
#           ),
#         width = 0.2, 
#         position = position_dodge(0.7)
#       ) +
#       common_y +
#       facet_grid(
#         rows = vars(Metric),
#         col = vars(Ancestry),
#         labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                             Metric = label_value)
#       ) +
#       labs(
#           # title = paste(gsub("_", " ", comparison)),
#           x = "Prediction",
#           y = "Correlation coefficient"
#       ) +
#       geom_text(
#         aes(
#           x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
#           y = Inf,  # Position the text at the top of the plot
#           label = paste("Perm. test,", "p.adj =", format(p_adjusted, digits = 3))
#           ),
#         size = 3,    
#         vjust = 1.5,   # Align text to the top (since we're using Inf for y)
#         hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
#         inherit.aes = FALSE  
#     ) 

#     # Add the plot to the list with a unique name
#     ggplot_list_dge[[paste0("prop_", prop)]] <- prop_plot
#   }
#   # Combine with patchwork 
#   combined_dge_bar_plot <- wrap_plots(ggplot_list_dge, ncol = 2)
#   # Combine ancestry wise
#   all_combined_plots[[ancestry]] <- combined_dge_bar_plot
# }

# final_patchwork_plot <- wrap_plots(all_combined_plots, ncol = 3) + 
#   plot_annotation(title = "DGE Robustness across ancestries") & 
#   theme(plot.title = element_text(hjust = 0.5))

# ggsave(
#   filename = "Plot_dge_robustness.pdf", 
#   plot = final_patchwork_plot, 
#   path = path_to_save_location, 
#   width = 15, height = 7.3
# )





# all_combined_plots <- list()
# for (ancestry in ancestries){
#   ggplot_list_dge = list()
#   for (prop in proportions) {
#     # Filter by proportion
#     filtered_dge <- summarized_metric_ml |> filter(
#       Proportion == prop,
#       Ancestry == ancestry
#       )

#     # Create lable to showcase number of samples
#     n_inf_label <- filtered_dge |>
#       ungroup() |>  
#       distinct(Ancestry, n_inf_ancestry) |>
#       mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
#       select(Ancestry, label) |>
#       deframe()
    
#     # Correlation axis from 0 to 1
#     common_y <- scale_y_continuous(
#       limits = c(0, 1.2), 
#       breaks = c(0, 0.5, 1))
    
#     # Create plot
#     prop_plot <- filtered_dge |>
#       ggplot(
#         aes(
#           x = fct_rev(Prediction),
#           y = mean_value
#           )
#       ) +
#       geom_bar(
#         stat = "identity", 
#         width = 0.7
#       ) +
#       geom_errorbar(
#         aes(
#           ymin = mean_value - sd_value, 
#           ymax = mean_value + sd_value
#           ),
#         width = 0.2, 
#         position = position_dodge(0.7)
#       ) +
#       common_y +
#       facet_grid(
#         rows = vars(Algorithm),
#         col = vars(Ancestry),
#         labeller = labeller(Ancestry = as_labeller(n_inf_label), 
#                             Algorithm = label_value)
#       ) +
#       labs(
#           # title = paste(gsub("_", " ", comparison)),
#           x = "Prediction",
#           y = "Correlation coefficient"
#       ) +
#       geom_text(
#         aes(
#           x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
#           y = Inf,  # Position the text at the top of the plot
#           label = paste("Perm. test,", "p.adj =", format(p_adjusted, digits = 3))
#           ),
#         size = 3,    
#         vjust = 1.5,   # Align text to the top (since we're using Inf for y)
#         hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
#         inherit.aes = FALSE  
#     ) 

#     # Add the plot to the list with a unique name
#     ggplot_list_dge[[paste0("prop_", prop)]] <- prop_plot
#   }
#   # Combine with patchwork 
#   combined_dge_bar_plot <- wrap_plots(ggplot_list_dge, ncol = 2)

#   # Combine ancestry wise
#   all_combined_plots[[ancestry]] <- combined_dge_bar_plot
# }

# final_patchwork_plot <- wrap_plots(all_combined_plots, ncol = 3) + 
#   plot_annotation(title = "ML Robustness across ancestries") & 
#   theme(plot.title = element_text(hjust = 0.5))

# ggsave(
#   filename = "Plot_ml_robustness.pdf", 
#   plot = final_patchwork_plot, 
#   path = path_to_save_location, 
#   width = 15, height = 7.3
# )






# # --------------------------------------------------------------------------------------------
# # Analysis: Interactions
# analysis_suffix <- "interactions"
# match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# # Extract files
# all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
# match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)
# # Print the match folders
# print("Matched folders:")
# print(match_vscratch_dir)

# # Check if there were matching folders
# if (length(match_vscratch_dir) == 0) {
#   message("No matching folders found.")
# }

# # Load interactions (Volcano plot)
# interactions <- data.frame()
# for (folder in match_vscratch_dir){
#     interaction_file <- file.path(folder, "Interactions.csv")
  
#     # Check if files exists
#     if (file.exists(interaction_file)) {
#       # Read data
#       interaction_data <- fread(interaction_file)
#       # Append data
#       interactions <- bind_rows(interactions, interaction_data)
#     } else {
#       warning("File does not exist: ", interaction_file)
#     }
# }

# # Volcano plot
# interactions <- filter(interactions, str_detect(coef, "\\."))
# logFC_threshold <- 1
# interactions <-  interactions |> 
#   mutate(Significant = adj.P.Val < 0.05 & abs(logFC) > logFC_threshold) |>
#   mutate(Ancestry = toupper(sub("\\..*", "", coef)))

# volcano_plot <- interactions |>
#   ggplot(
#     aes(
#       x = logFC,
#       y = -log10(adj.P.Val),
#       color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
#     )
#   ) +
#   geom_point() +
#   facet_grid(cols = vars(Ancestry)) + 
#   geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
#   scale_color_manual(
#     values = c("TRUE" = "red", "FALSE" = "grey"),
#   ) +
#   labs(
#     title = "Significant interactions across ancestries"
#   ) +
#   # geom_text_repel(
#   #   data = subset(volcano_data, Significant),  
#   #   aes(label = Feature),                     
#   #   size = 3,                                   
#   #   color = "black",
#   #   # max.overlaps = Inf                              
#   # ) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#     )

# # Save the plot
# ggsave(filename = "Volcano_plot.pdf", 
#        plot = volcano_plot, 
#        path = path_to_save_location, 
#        width = 7, height = 3)


# # Load enrichment results
# enrichment <- data.frame()
# for (folder in match_vscratch_dir){
#     interaction_file <- file.path(folder, "Interactions_FGSEA_enrichment.csv")
  
#     # Check if files exists
#     if (file.exists(interaction_file)) {
#       # Read data
#       interaction_data <- fread(interaction_file)
#       # Append data
#       enrichment <- bind_rows(enrichment, interaction_data)
#     } else {
#       warning("File does not exist: ", interaction_file)
#     }
# }

# FGSEA_plot_1 <- enrichment |>
#   filter(
#     Top30 == TRUE,
#     Database == "MSigDB Hallmark 2020") |>
#   ggplot(
#     aes(
#       x = Ancestry,
#       y = pathway,
#       color = NES,
#       size = pmin(-log10(padj), 5)
#     )
#   ) +
#   geom_point() +
#   scale_size_binned(
#     range = c(1, 3)    
#   ) +
#   scale_color_gradient2(
#     high = "red", 
#     mid = "white", 
#     low = "blue"
#   ) +
#   labs(
#     title = unique(interactions$Database[interactions$Database == "MSigDB Hallmark 2020"]),
#     y = "Pathway",
#     size = "P.Val.Adj."
#   )

# # Save the plot
# ggsave(filename = "Interactions_1.pdf", 
#        plot = FGSEA_plot_1, 
#        path = path_to_save_location, 
#        width = 5, height = 6)


# FGSEA_plot_2 <- interactions |>
#   filter(
#     Top10 == TRUE,
#     Database != "MSigDB Hallmark 2020") |>
#   ggplot(
#     aes(
#       x = Ancestry,
#       y = pathway,
#       color = NES,
#       size = pmin(-log10(padj), 5)
#     )
#   ) +
#   geom_point() +
#   scale_size_binned(
#     range = c(1, 3)    
#   ) +
#   scale_color_gradient2(
#     high = "red", 
#     mid = "white", 
#     low = "blue"
#   ) +
#   labs(
#     title = unique(interactions$Database[interactions$Database != "MSigDB Hallmark 2020"]),
#     y = "Pathway",
#     size = "P.Val.Adj."
#   )

# # Save 
# ggsave(filename = "Interactions_2.pdf", 
#        plot = FGSEA_plot_2, 
#        path = path_to_save_location, 
#        width = 8, height = 4)



# # Bigger Patchwork plot
# generalizability_plot <- dge_bar_ml_bar +
#   plot_annotation(title = "Generalizability across ancestries")

# patchwork1 <- (generalizability_plot | (performance_plot / volcano_plot)) +
#   plot_layout(widths = c(2, 1)) 

# # Save 
# ggsave(filename = "Patchwork_1.pdf", 
#        plot = patchwork1, 
#        path = path_to_save_location, 
#        width = 15, height = 7.5)
