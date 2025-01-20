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
print("Evoking R.")

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
ml_metric <- ml_metric |>  
    pivot_longer(
        cols = c(LogisticRegression, RandomForestClassifier),
        values_to = 'Value',
        names_to = 'Algorithm'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))


dge_metric <- dge_metric |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))


# Define common scale for y-axis
common_scale_y <- scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
)

# Create ancestry labels with n_ancestry
ancestry_labels <- ml_metric |>
  distinct(Ancestry, n_inf_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()  # Converts to a named vector

# Visualize:
ml_metric_plot <- ml_metric |>
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
        rows = vars(Algorithm),
        labeller = labeller(
            Ancestry = as_labeller(ancestry_labels),
            Metric = label_value
        )
    ) +
    labs(
        x = "Prediction",
        y = "ROC AUC"
    )

# Save
# Save the image
ggsave(filename = "Plot_ml.pdf", 
       plot = ml_metric_plot, 
       path = setup$output_directory, 
       width = 5, height = 5)



dge_metric_plot <- dge_metric |> 
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
    labs(
        x = "Prediction",
        y = "Correlation coefficient"
    )

# Save
ggsave(filename = "Plot_dge.pdf", 
       plot = dge_metric_plot, 
       path = setup$output_directory, 
       width = 5, height = 5)



