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
ml_metric <- fread(file.path(setup$output_directory, "Contrast_metric_ml.csv")) 
dge_metric <- fread(file.path(setup$output_directory, "Contrast_metric_dge.csv")) 


# Change format (pivot longer)
ml_metric <- ml_metric |>  
    pivot_longer(
        cols = c(ROC_AUC, Accuracy),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    mutate(Prediction = fct_rev(Prediction))


correlation_dge_metric <- dge_metric |> 
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

# Create threshold label


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
    rows = vars(Algorithm, Metric),
    labeller = labeller(
      Ancestry = as_labeller(ancestry_labels)
      )
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  )

# Save
# Save the image
ggsave(filename = "Plot_phenotype_prediction.pdf", 
       plot = ml_metric_plot, 
       path = setup$output_directory, 
       width = 5, height = 5)



dge_correlation_plot <- correlation_dge_metric |> 
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
ggsave(filename = "Plot_logFC_correlation.pdf", 
       plot = dge_correlation_plot, 
       path = setup$output_directory, 
       width = 5, height = 5)


# # DGE prediction
# prediction_dge_metric <- dge_metric |>
#     pivot_longer(
#         cols = c(RMSE, R2),
#         values_to = 'Value',
#         names_to = 'Metric'
#     ) |> 
#     mutate(Prediction = fct_rev(Prediction))

# # Visualize
# dge_prediction_plot <- prediction_dge_metric |>
#     ggplot(
#         aes(
#             x = Prediction,
#             y = Value
#         )
#     ) +
#     geom_col() +
#     common_scale_y +
#     facet_grid(
#         cols = vars(Ancestry),
#         rows = vars(Metric),
#         labeller = labeller(
#             Ancestry = as_labeller(ancestry_labels),
#             Metric = label_value
#         )
#     ) +
#     labs(
#         x = "Prediction",
#         y = "Y"
#     )

# ggsave(filename = "Plot_logFC_prediction.pdf", 
#        plot = dge_prediction_plot, 
#        path = setup$output_directory, 
#        width = 5, height = 5)

