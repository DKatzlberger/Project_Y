suppressPackageStartupMessages({
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    # Visualization
    library(patchwork)
})

# Custom functions
source("r_utils.R")

# Here starts the script
print("Envoking R.")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  YAML_FILE <- args[1]  
  # Check if it's a valid YAML file
  is_yml_file(YAML_FILE)
  setup <- yaml.load_file(YAML_FILE)

} else {
  # Dev settings if no command-line argument provided
  YAML_FILE <- "dev_settings.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(YAML_FILE)
}


# Here starts the script
print('Start visualization.')

# Load metric dataframes
ml_metric <- fread(file.path(setup$output_directory, 'Metric_ml.csv')) 
dge_metric <- fread(file.path(setup$output_directory, 'Metric_dge.csv')) 


# Change format (pivot longer)
metric <- ml_metric |>  
    pivot_longer(
        cols = c(Accuracy, F1, ROC_AUC),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))

# Plots
# Define common scale for y-axis
common_scale_y <- scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
)

# Create ancestry labels with n_ancestry
ancestry_labels <- metric |>
  distinct(Ancestry, n_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()  # Converts to a named vector

# Create threshold labels
threshold_labels <- metric |>
  distinct(Threshold) |>
  mutate(label = paste0("Thr: ", Threshold)) |>
  select(Threshold, label) |>
  deframe()  # Converts to a named vector

# Plot F1
f1_plot <- metric |> 
    filter(Metric == 'F1') |> 
    ggplot(
        aes(
            x = Prediction,
            y = Value
        )
    ) +
    geom_col() +
    common_scale_y +
    facet_grid(
        cols = vars(Ancestry),
        rows = vars(Threshold),
        labeller = labeller(
            Ancestry = as_labeller(ancestry_labels),
            Threshold = as_labeller(threshold_labels)
            )
    ) +
    ggtitle('F1') 

# Plot Accuracy
acc_plot <- metric |> 
    filter(Metric == 'Accuracy') |> 
    ggplot(
        aes(
            x = Prediction,
            y = Value
        )
    ) +
    geom_col() +
    common_scale_y +
    facet_grid(
        cols = vars(Ancestry),
        rows = vars(Threshold),
        labeller = labeller(
            Ancestry = as_labeller(ancestry_labels),
            Threshold = as_labeller(threshold_labels)
            )
    ) +
    ggtitle('Accuracy')

# Plot ROC AUC
auc_plot <- metric |> 
    filter(Metric == 'ROC_AUC') |> 
    select(!Threshold) |> 
    unique() |> 
    ggplot(
        aes(
            x = Prediction,
            y = Value
        )
    ) +
    geom_col() +
    common_scale_y +
    facet_grid(
        cols = vars(Ancestry),
        labeller = labeller(
            Ancestry = as_labeller(ancestry_labels)
        )
    ) +
    ggtitle('ROC AUC') 
    
# Patchwork 
combined <- f1_plot + acc_plot + auc_plot +
  plot_layout(ncol = 2)

# Save
ggsave(file.path(setup$output_directory, 'Generalizability_ml.pdf'), height = 10, width = 10)

# Plot correlation of logFCs
metric <- dge_metric |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))

log_correlation <- metric |> 
    ggplot(
        aes(
            x = Prediction,
            y = Value
        )
    ) +
    geom_col() +
    common_scale_y +
    facet_grid(
        cols = vars(Ancestry),
        rows = vars(Metric),
        labeller = labeller(
            Ancestry = as_labeller(ancestry_labels),
            Metric = label_value
        )
    ) +
    ggtitle('Correlation') 

# Save
ggsave(file.path(setup$output_directory, 'Generalizability_dge.pdf'), height = 10, width = 10)


