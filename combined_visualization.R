# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    }
)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  YAML_FILE <- args[1]  # First argument is the YAML file path
} else {
  # Dev settings
  YAML_FILE <- "job_settings.yml"
  print("Running interactive mode for development.\n")
}

# Load setup file
setup <- yaml.load_file(YAML_FILE)

# Location where I store my files
vscratch_dir <- "data/runs"

# Tag is used to know which data is used
tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
train_ancestry <- toupper(setup$classification$train_ancestry)
infer_ancestry <- toupper(setup$classification$infer_ancestry)

# Combine
modified_path <- paste0(comparison, "_", train_ancestry, "_to_", infer_ancestry)

# Path to the directories (without seed information)
# Used for the pattern
directory <- file.path(vscratch_dir, modified_path)

# Extract all folders matching the pattern in the vscratch directory
all_folders <- list.dirs(vscratch_dir, full.names = TRUE, recursive = FALSE)
matching_folders <- grep(modified_path, all_folders, value = TRUE)

# Combining metric csv files 
metric_dge <- data.frame()
metric_ml <- data.frame()
for (folder in matching_folders){
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}

# Visualize all seeds (DGE and ML)
# Stats
# Prepare data format
summarized_dge <- metric_dge |> 
    mutate(Prediction = fct_rev(Prediction)) |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    group_by(Ancestry, Status, Prediction, Metric) |>  
    summarize(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n())
  ) 

# ML
# Prepare data format
summarized_ml <- metric_ml |> 
    # select(c(Status, Ancestry, Seed, Prediction, ROC_AUC)) |> 
    mutate(Prediction = fct_rev(Prediction)) |> 
    pivot_longer(
        cols = c(ROC_AUC),
        values_to = 'Value',
        names_to = 'Metric'
    ) |> 
    group_by(Ancestry, Status, Prediction, Metric) |>  
    summarize(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n())
  ) 

# Space for some shared variables (plotting)

# Plot of correlation 
DGE_plot <- summarized_dge |> 
    ggplot(
        aes(
            x = Prediction,
            y = mean_value
        )
    ) +
    geom_bar(
        stat = "identity", 
        width = 0.7
    ) +
    geom_errorbar(
        aes(
            ymin = mean_value - se_value, 
            ymax = mean_value + se_value
            ), 
        width = 0.2, position = position_dodge(0.7)
    ) +
    facet_grid(
        rows = vars(Metric),
        col = vars(Ancestry)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "Correlation coefficient"
    )

# Plot of ml
ML_plot <- summarized_ml |> 
        ggplot(
        aes(
            x = Prediction,
            y = mean_value
        )
    ) +
    geom_bar(
        stat = "identity", 
        width = 0.7
    ) +
    geom_errorbar(
        aes(
            ymin = mean_value - se_value, 
            ymax = mean_value + se_value
            ), 
        width = 0.2, position = position_dodge(0.7)
    ) +
    facet_grid(
        # rows = vars(Metric),
        col = vars(Ancestry)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "ROC AUC"
    )

# Save the data (create directory)
# Make a combined directory for the comparison (based on the var: comparison)
base_dir <- 'data/combined_runs'
combined_dir <- file.path(base_dir, comparison)
# Create directory if it does not exists
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

save_directory <- file.path(combined_dir, modified_path)
# Create directory if it does not exists
if (!dir.exists(save_directory)) {
  dir.create(save_directory)
}

# Save files and plots
fwrite(metric_dge, file.path(save_directory, 'Metric_dge.csv'))
fwrite(metric_ml, file.path(save_directory, 'Metric_ml.csv'))

# Save the image
ggsave(filename = 'Plot_dge.pdf', plot = DGE_plot, 
       path = save_directory, width = 8, height = 6
       )

ggsave(filename = 'Plot_ml.pdf', plot = ML_plot, 
       path = save_directory, width = 8, height = 6
       )

# Save the R object
saveRDS(DGE_plot, file = file.path(save_directory, 'Plot_dge.rds'))
saveRDS(ML_plot, file = file.path(save_directory, 'Plot_ml.rds'))

# Combine model weights
# Read the files 
combined_weights <- data.frame()
for (folder in matching_folders){
    weights_file <- file.path(folder, "Weights.csv")

    # Load and append DGE data
    weights <- fread(weights_file) 
    combined_weights <- bind_rows(combined_weights, weights) 
}
# # Take the mean (summarize)
# summarized_weights = as_tibble(colMeans(combined_weights), rownames = 'Feature')

# # Function to calculate SEM
# calculate_sem <- function(x) {
#   return(sd(x) / sqrt(length(x)))  # Standard Error of the Mean (SEM)
# }

# # Calculate SEM for each column and convert to a tibble
# sem_weights <- apply(combined_weights, 2, calculate_sem)

# # Add SEM to summarized_weights
# summarized_weights$se_value <- sem_weights

# Save
fwrite(combined_weights, file.path(save_directory, 'Weights.csv'))






















