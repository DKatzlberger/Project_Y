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
    }
)
# Source custom functions
source("r_utils.R")

# Combine all ancestries for a comparison
vscratch_dir_in = file.path("data", "combined_runs")
# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders

# Create pattern
comparison <- "PanCanAtlas_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma"
analysis_suffix <- "cross_ancestry"
match_pattern <- paste0(comparison, ".*", analysis_suffix, "$")

# Extract files
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)
# Print the match folders
print("Matched folders:")
print(match_vscratch_dir)

# Check if there were matching folders
if (length(match_vscratch_dir) == 0) {
  message("No matching folders found.")
}

# Save the results of the analysis
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- file.path("data", "combined_ancestries")
path_to_save_location <- file.path(vscratch_dir_out, comparison)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Combining metric csv files 
# Each folder is one seed
metric_dge <- data.frame()
metric_ml <- data.frame()
for (folder in match_vscratch_dir){
    dge_file <- file.path(folder, "Summarized_metric_dge.csv")
    ml_file <- file.path(folder, "Summarized_metric_ml.csv")

    # Load and append DGE data for each seed
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data for each seed
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}

# Pvalue correction
# dge
unique_p_values_dge <- metric_dge |> pull(p_value) |> unique()
# Benjamini-Hochberg correction
p_adjusted_bh_dge <- p.adjust(unique_p_values_dge, method = "BH")
# Create a named vector to map the adjusted p-values back to the data frame
p_adjusted_named_dge <- setNames(p_adjusted_bh_dge, unique_p_values_dge)
# Add 
metric_dge <- metric_dge |>
   mutate(p_adjusted = p_adjusted_named_dge[as.character(p_value)])

# ml
unique_p_values_ml <- metric_ml |> pull(p_value) |> unique()
# Benjamini-Hochberg correction
p_adjusted_bh_ml <- p.adjust(unique_p_values_ml, method = "BH")
# Create a named vector to map the adjusted p-values back to the data frame
p_adjusted_named_ml <- setNames(p_adjusted_bh_ml, unique_p_values_ml)
# Add 
metric_ml <- metric_ml |>
   mutate(p_adjusted = p_adjusted_named_ml[as.character(p_value)])

# Save the metric data frames
fwrite(metric_dge, file.path(path_to_save_location, "Metric_dge.csv"))
fwrite(metric_ml, file.path(path_to_save_location, "Metric_ml.csv"))

# Visualization:
ancestries <- unique(metric_dge$Ancestry)
# dge
ggplot_list_dge = list()
for (ancestry in ancestries){
  filtered_metric <- filter(metric_dge, Ancestry == ancestry)
  # Create lable to showcase number of samples
  n_inf_label <- filtered_metric |>
    ungroup() |>  
    distinct(Ancestry, n_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_ancestry, ")")) |>
    select(Ancestry, label) |>
    deframe()

  # Correlation axis from 0 to 1
  common_y <- scale_y_continuous(
    limits = c(0, 1.2), 
    breaks = c(0, 0.5, 1))

  plot <- filtered_metric |>
    ggplot(
      aes(
        x = fct_rev(Prediction),
        y = mean_value
        )
    ) +
    geom_bar(
      stat = "identity", 
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        ymin = mean_value - sd_value, 
        ymax = mean_value + sd_value
        ),
      width = 0.2, 
      position = position_dodge(0.7)
    ) +
    common_y +
    facet_grid(
      rows = vars(Metric),
      col = vars(Ancestry),
      labeller = labeller(Ancestry = as_labeller(n_inf_label), 
                          Metric = label_value)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "Correlation coefficient"
    ) +
    geom_text(
      aes(
        x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
        y = Inf,  # Position the text at the top of the plot
        label = paste("Perm. test,", "p =", format(p_adjusted, digits = 3))
        ),
      size = 4,    
      vjust = 2,   # Align text to the top (since we're using Inf for y)
      hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
      inherit.aes = FALSE  
  ) 

  # Add the plot to the list with a unique name
  ggplot_list_dge[[paste0("ancestry_", ancestry)]] <- plot
}

# Combine with patchwork 
combined_dge_plot <- wrap_plots(ggplot_list_dge, ncol = 2)
# Save the image
ggsave(filename = "Plot_dge.pdf", 
       plot = combined_dge_plot, 
       path = path_to_save_location, 
       width = 10, height = 10)

# ml
# ml 
ggplot_list_ml = list()
for (ancestry in ancestries) {
  filtered_metric <- filter(metric_ml, Ancestry == ancestry)

  # Create lable to showcase number of samples
  n_inf_label <- filtered_metric |>
    ungroup() |>  
    distinct(Ancestry, n_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_ancestry, ")")) |>
    select(Ancestry, label) |>
    deframe()
  
  # Correlation axis from 0 to 1
  common_y <- scale_y_continuous(
    limits = c(0, 1.2), 
    breaks = c(0, 0.5, 1))
  
  # Create plot
  plot <- filtered_metric |>
    ggplot(
      aes(
        x = fct_rev(Prediction),
        y = mean_value
        )
    ) +
    geom_bar(
      stat = "identity", 
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        ymin = mean_value - sd_value, 
        ymax = mean_value + sd_value
        ),
      width = 0.2, 
      position = position_dodge(0.7)
    ) +
    common_y +
    facet_grid(
      # rows = vars(Metric),
      col = vars(Ancestry),
      labeller = labeller(Ancestry = as_labeller(n_inf_label), 
                          Metric = label_value)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "ROC AUC"
    ) +
    geom_text(
      aes(
        x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
        y = Inf,  # Position the text at the top of the plot
        label = paste("Perm. test,", "p =", format(p_adjusted, digits = 3))
        ),
      size = 4,    
      vjust = 2,   # Align text to the top (since we're using Inf for y)
      hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
      inherit.aes = FALSE  
  ) 

  # Add the plot to the list with a unique name
  ggplot_list_ml[[paste0("ancestry_", ancestry)]] <- plot
}

# Combine with patchwork 
combined_ml_plot <- wrap_plots(ggplot_list_ml, ncol = 2)
# Save the image
ggsave(filename = "Plot_ml.pdf", 
       plot = combined_ml_plot, 
       path = path_to_save_location, 
       width = 10, height = 10)