suppressPackageStartupMessages({
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    # Visualization
    library(patchwork)
    library(httpgd)
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv('/opt/conda/envs/ancestry/bin/python')
})

# Custom functions
source('r_utils.R')

# Here starts the script
print('Envoking R.')
# Load the settings (They are a command line input)
args = commandArgs(trailingOnly = TRUE)
# Check that it is a yaml file and exists
is_yml_file(args[1])
setup <- yaml.load_file(args[1])

# Run file for binary metric
args  <- args[1]
python_script_path  <- file.path(getwd(), 'metric_calculation.py')

py_run_string(sprintf("import sys; sys.argv = ['%s', '%s']", 
                      python_script_path, paste(args, collapse="', '")))

py_run_file(python_script_path)

# Here starts the script
print('Start visualization.')

# Change format (pivot longer)
metric <- py$r_metric |>  
    pivot_longer(
        cols = c(Accuracy, F1, `ROC AUC`),
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
    filter(Metric == 'ROC AUC') |> 
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
f1_plot + acc_plot + auc_plot +
  plot_layout(ncol = 2)

# Save
ggsave(file.path(setup$output_directory, 'Validation_ml.pdf'), height = 12, width = 12)


