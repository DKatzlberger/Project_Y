suppressPackageStartupMessages({
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    # Visualization
    library(patchwork)
})

# Custom functions
source('r_utils.R')

# Here starts the script
print('Envoking R.')

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
test_metric <- fread(file.path(setup$output_directory, 'Metric_test.csv')) |> 
    mutate(Prediction = 'Subset')
inf_metric  <- fread(file.path(setup$output_directory, 'Metric_inf.csv')) |> 
    mutate(Prediction = 'Ancestry')

# Combine for visulization
metric <- 
    bind_rows(test_metric, inf_metric)


# Change format (pivot longer)
metric <- metric |>  
    pivot_longer(
        cols = c(Accuracy, F1, ROC_AUC),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))

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
    facet_grid(
        cols = vars(Ancestry),
        rows = vars(Threshold)
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
    facet_grid(
        cols = vars(Ancestry),
        rows = vars(Threshold)
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
    facet_grid(
        cols = vars(Ancestry),
    ) +
    ggtitle('ROC AUC') 
    
# Patchwork 
combined <- f1_plot + acc_plot + auc_plot +
  plot_layout(ncol = 2)

# Save
ggsave(file.path(setup$output_directory, 'Validation_ml.pdf'), height = 10, width = 10)


