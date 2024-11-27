# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Statistics 
    library(coin)
    }
)
# Source custom functions
source("r_utils.R")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  yaml_file <- args[1]  # First argument is the YAML file path
} else {
  # Dev settings
  yaml_file <- "job_settings.yml"
  print("Running interactive mode for development.\n")
}

# Load setup file
setup <- yaml.load_file(yaml_file)

# Location where I store my files
vscratch_dir <- "data/runs"

# Tag is used to know which data is used
tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
train_ancestry <- toupper(setup$classification$train_ancestry)
infer_ancestry <- toupper(setup$classification$infer_ancestry)

# Combine (modified path is used for the pattern)
modified_path <- paste0(comparison, "_", train_ancestry, "_to_", infer_ancestry)

# Is used when there is a dir for all runs otherwise it is vscratch_dir
# Directory to run folders for comparison
directory <- file.path(vscratch_dir, comparison)

# Extract all folders matching the pattern in the vscratch directory
all_folders <- list.dirs(directory, full.names = TRUE, recursive = FALSE)
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
# Add information to the dataframe
metric_dge <- metric_dge |> mutate(Comparison = comparison)
metric_ml <- metric_ml |> mutate(Comparison = comparison)

# Outlier test
# With hampel filter 3 * MAD confidence interval
hampel_filter <- function(x, t0 = 3) {
  # x: Numeric vector
  # t0: Threshold (default is 3 standard deviations equivalent)
  med <- median(x, na.rm = TRUE)  # Median of the data
  mad_val <- mad(x, na.rm = TRUE)  # MAD of the data
  threshold <- t0 * mad_val  # Threshold based on MAD and t0
  
  # Identify outliers
  abs(x - med) > threshold
}

# Add outlier column for correlation values
metric_dge <- metric_dge |>
  group_by(Status) |>
  mutate(
    Pearson_outlier = hampel_filter(Pearson),  # Hampel filter for Pearson
    Spearman_outlier = hampel_filter(Spearman),  # Hampel filter for Spearman
  )

# Add otlier column to check if ROC_AUC values contain outliers
metric_ml <- metric_ml |>
  group_by(Status) |>
  mutate(
    ROC_AUC_outlier = hampel_filter(ROC_AUC),  # Hampel filter for ROC AUC
  )

# Differential gene expression analysis:

# Test difference in distribution of correlated logFC  
# 1. with Kolmogorow-Smirnow
# Non-parametric
# H0: The samples come from the same probability distribution
# H1: The samples are from different probability distribution
ks_pearson <- ks_test(data = metric_dge, group_col = "Status", value_col = "Pearson")
ks_spearman <- ks_test(data = metric_dge, group_col = "Status", value_col = "Spearman")

# 2. with permutation testing (returns a p-value of 1 if there is no variance between groups)
# Non-parametric
# Without assumption
p_pearson <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Pearson")
p_spearman <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Spearman")

# Summarize data (taking mean of correlation values)
summarized_dge <- metric_dge |> 
    mutate(Prediction = fct_rev(Prediction)) |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = "Value",
        names_to = "Metric"
    ) |> 
    group_by(Ancestry, Status, Prediction, Metric, n_ancestry) |>  
    summarize(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n())
  ) |>
  mutate(p_value = ifelse(Metric == "Pearson", p_pearson, p_spearman))

# Machine learning:

# Test difference in distribution of correlated logFC  
# 1. with Kolmogorow-Smirnow
# Non-parametric
# H0: The samples come from the same probability distribution
# H1: The samples are from different probability distribution
ks_roc <- ks_test(data = metric_ml, group_col = "Status", value_col = "ROC_AUC")

# 2. with permutation testing
# Non-parametric
# Without assumption
p_roc <- permutation_test(data = metric_ml, group_col = "Status", value_col = "ROC_AUC")

# Prepare data format
summarized_ml <- metric_ml |> 
    # select(c(Status, Ancestry, Seed, Prediction, ROC_AUC)) |> 
    mutate(Prediction = fct_rev(Prediction)) |> 
    pivot_longer(
        cols = c(ROC_AUC),
        values_to = "Value",
        names_to = "Metric"
    ) |> 
    group_by(Ancestry, Status, Prediction, Metric, n_ancestry) |>  
    summarize(
    mean_value = mean(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n())
  ) |>
  mutate(p_value = p_roc)

# Plot differential gene expression analysis and machine learning
# Space for some shared variables (plotting)

# Labeller to annotate number of samples per ancestry
ancestry_labels <- summarized_dge |>
  ungroup() |>  
  distinct(Ancestry, n_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()

# Correlation axis from 0 to 1
common_y <- scale_y_continuous(
    limits = c(0, 1.1), 
    breaks = c(0, 0.5, 1))

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
    common_y +
    facet_grid(
        rows = vars(Metric),
        col = vars(Ancestry),
        labeller = labeller(Ancestry = as_labeller(ancestry_labels), 
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
      label = paste("Perm. test,", "p =", format(p_value, digits = 3)),  # Add text for the p-value source
    ),
    size = 4,    # Adjust text size
    vjust = 2,   # Align text to the top (since we're using Inf for y)
    hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
    inherit.aes = FALSE  # Don't inherit the default aesthetics
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
    common_y +
    facet_grid(
        # rows = vars(Metric),
        col = vars(Ancestry),
        labeller = labeller(Ancestry = as_labeller(ancestry_labels))
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
      label = paste("Perm. test,", "p =", format(p_value, digits = 3)),  # Add text for the p-value source
    ),
    size = 4,    # Adjust text size
    vjust = 2,   # Align text to the top (since we're using Inf for y)
    hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
    inherit.aes = FALSE  # Don't inherit the default aesthetics
  )


# Save the data (create directory)
# Make a combined directory for the comparison (based on the var: comparison)
base_dir <- "data/combined_runs"
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
fwrite(metric_dge, file.path(save_directory, "Metric_dge.csv"))
fwrite(metric_ml, file.path(save_directory, "Metric_ml.csv"))

# Save the image
ggsave(filename = "Plot_dge.pdf", plot = DGE_plot, 
       path = save_directory, width = 5, height = 5
       )

ggsave(filename = "Plot_ml.pdf", plot = ML_plot, 
       path = save_directory, width = 5, height = 5
       )

# Save the R object
saveRDS(DGE_plot, file = file.path(save_directory, "Plot_dge.rds"))
saveRDS(ML_plot, file = file.path(save_directory, "Plot_ml.rds"))

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
fwrite(combined_weights, file.path(save_directory, "Weights.csv"))

# TODO - Create overall directory in runs to move single run files






















