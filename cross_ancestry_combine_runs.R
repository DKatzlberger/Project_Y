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
  yaml_file <- args[1] 
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
  # Vscratch directory
  vscratch_dir_in = file.path("data", "runs")
  # When the script is run from the command line then 'output_directory' is given
  # The pattern to extract all matchin directories is extracted from 'output_directory'
  output_path = setup$output_directory
  match_pattern <- sub("_\\d+$", "", sub(".*/", "", output_path))
} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "PanCanAtlas_RSEM_basal_vs_non-basal_EUR_to_ADMIX.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

  # Construction:
  # Vscratch_dir is the place where the files are stored
  vscratch_dir_in = file.path("data", "runs")
  # Tag is used to specify which data it is e.g. TCGA, NIAGADS
  tag <- setup$tag
  # Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
  comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
  # Analysis_name specifies what was analysed
  # E.g. comparsion of ancestries EUR_to_AMR, subsetting_EUR
  # This often is modified depending which analysis
  train_ancestry <- toupper(setup$classification$train_ancestry)
  infer_ancestry <- toupper(setup$classification$infer_ancestry)
  analysis_name  <- paste0(train_ancestry, "_to_", infer_ancestry, "_cross_ancestry")
  # Combination of components to create the 'match_pattern'
  # The 'match_pattern' is used as pattern to extract all folders in the vscratch dir
  match_pattern <- paste0(comparison, "_", analysis_name)
}

# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders
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
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
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
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data for each seed
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data for each seed
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}
# Save
fwrite(metric_dge, file.path(path_to_save_location, "Metric_dge.csv"))
fwrite(metric_ml, file.path(path_to_save_location, "Metric_ml.csv"))

# Similarity of obervations across seeds
# Load observations
subset_obs <- data.frame()
for (folder in match_vscratch_dir){

    obs_test_file <- file.path(folder, "Obs_test.yml")

    # Load and append obs for each seed
    obs_test_data <- yaml.load_file(obs_test_file)
    # Seed
    seed <- sub(".*_(\\d+)$", "\\1", folder)
    # Make dataframe
    obs_test_df <- data.frame(Samples = obs_test_data,
                              Seed = seed)
    # Combine
    subset_obs <- bind_rows(subset_obs, obs_test_df)                   
}

# Group samples by seed
seeds_samples <- subset_obs |>
  group_by(Seed) |>
  summarise(Samples = list(Samples), .groups = "drop")

# Compute pairwise Jaccard indices and overlap coefficients
pairwise_metrics <- combn(seeds_samples$Samples, 2, function(pair) {
  intersection <- length(intersect(pair[[1]], pair[[2]]))
  union <- length(union(pair[[1]], pair[[2]]))
  min_size <- min(length(pair[[1]]), length(pair[[2]]))
  
  jaccard <- intersection / union
  overlap <- intersection / min_size
  
  c(jaccard = jaccard, overlap = overlap)
}, simplify = FALSE)

# Convert pairwise metrics to a data frame
metrics_df <- do.call(rbind, pairwise_metrics) |>
  as.data.frame() |>
  setNames(c("Jaccard", "Overlap"))

# Compute average Jaccard Index and Overlap Coefficient
average_jaccard <- mean(metrics_df$Jaccard)
average_overlap <- mean(metrics_df$Overlap)

# Visualization
# Correlation axis from 0 to 1
common_y <- scale_y_continuous(
    limits = c(0, 1.1), 
    breaks = c(0, 0.5, 1))

# Create a dataframe for plotting
metrics_summary <- data.frame(
  Metric = c("Jaccard Index", "Overlap Coefficient"),
  Value = c(average_jaccard, average_overlap)
)

# Plot
sample_distribution_plot <- metrics_summary |>
  ggplot(
    aes(
      x = Metric,
      y = Value
    )
  ) +
  geom_col() +
  common_y + 
  labs(
    title = "Average pairwise overlap of subsets across seeds",
    x = "Metric",
    y = "Value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save
ggsave(filename = "Observation_overlap_subset.pdf", 
       plot = sample_distribution_plot, 
       path = path_to_save_location, 
       width = 5, height = 5
       )


# Add information to the dataframe
# comparison = paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
# metric_dge <- metric_dge |> mutate(Comparison = comparison)
# metric_ml <- metric_ml |> mutate(Comparison = comparison)

# Outlier test
# With hampel filter 3 * MAD confidence interval
# hampel_filter <- function(x, t0 = 3) {
#   # x: Numeric vector
#   # t0: Threshold (default is 3 standard deviations equivalent)
#   med <- median(x, na.rm = TRUE)  # Median of the data
#   mad_val <- mad(x, na.rm = TRUE)  # MAD of the data
#   threshold <- t0 * mad_val  # Threshold based on MAD and t0
  
#   # Identify outliers
#   abs(x - med) > threshold
# }

# # Add outlier column for correlation values
# metric_dge <- metric_dge |>
#   group_by(Status) |>
#   mutate(
#     Pearson_outlier = hampel_filter(Pearson),  # Hampel filter for Pearson
#     Spearman_outlier = hampel_filter(Spearman),  # Hampel filter for Spearman
#   )


# Add otlier column to check if ROC_AUC values contain outliers
# metric_ml <- metric_ml |>
#   group_by(Status) |>
#   mutate(
#     ROC_AUC_outlier = hampel_filter(ROC_AUC),  # Hampel filter for ROC AUC
#   )

# Differential gene expression analysis:

# Test difference in distribution of correlated logFC  
# 1. Kolmogorow-Smirnow
# Non-parametric
# H0: The samples come from the same probability distribution
# H1: The samples are from different probability distribution
# ks_pearson <- ks_test(data = metric_dge, group_col = "Status", value_col = "Pearson")
# ks_spearman <- ks_test(data = metric_dge, group_col = "Status", value_col = "Spearman")

# 2. with permutation testing (returns a p-value of 1 if there is no variance between groups)
# Non-parametric
# Without assumption
p_pearson <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Pearson")
p_spearman <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Spearman")

# Summarize data (taking mean of correlation values)
summarized_dge <- metric_dge |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = "Value",
        names_to = "Metric"
    ) |> 
    group_by(Ancestry, 
             Status, 
             Prediction, 
             Metric, 
             n_inf_ancestry
             ) |>  
    summarize(
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
  ) |>
  mutate(p_value = ifelse(Metric == "Pearson", p_pearson, p_spearman))

# Machine learning:
# 2. with permutation testing
# Non-parametric
# Without assumption
p_regression <- permutation_test(data = metric_ml, group_col = "Status", value_col = "LogisticRegression")
p_forest <- permutation_test(data = metric_ml, group_col = "Status", value_col = "RandomForestClassifier")

# Prepare data format
summarized_ml <- metric_ml |> 
    pivot_longer(
        cols = c(LogisticRegression, RandomForestClassifier),
        values_to = "Value",
        names_to = "Algorithm"
    ) |> 
    group_by(Ancestry, 
             Status, 
             Prediction, 
             Algorithm, 
             n_inf_ancestry
             ) |>  
    summarize(
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
  ) |>
  mutate(p_value = ifelse(Algorithm == "LogisticRegression", p_regression, p_forest))

# Save
fwrite(summarized_dge, file.path(path_to_save_location, "Summarized_metric_dge.csv"))
fwrite(summarized_ml, file.path(path_to_save_location, "Summarized_metric_ml.csv"))

# Visualization

# Labeller to annotate number of samples per ancestry
ancestry_labels <- summarized_dge |>
  ungroup() |>  
  distinct(Ancestry, n_inf_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()

# Axis
common_y <- scale_y_continuous(
  limits = c(0, 1.2), 
  breaks = c(0, 0.5, 1))

common_x <- scale_x_continuous(
  limits = c(-1.2, 1.2), 
  breaks = c(-1,-0.5,0, 0.5, 1))

# Bar plot dge
DGE_bar_plot <- summarized_dge |> 
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
      label = ifelse(
        is.na(p_value),  
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p =", format(p_value, digits = 3))  # Otherwise, display the p-value
      )
    ),
    size = 4,    # Adjust text size
    vjust = 2,   # Align text to the top (since we're using Inf for y)
    hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
    inherit.aes = FALSE  # Don't inherit the default aesthetics
  ) 

# Save the plot
ggsave(filename = "Plot_bar_dge.pdf", 
       plot = DGE_bar_plot, 
       path = path_to_save_location, 
       width = 5, height = 5
       )

# Density plot
DGE_density_plot <- metric_dge |>
  pivot_longer(
    cols = c(Pearson, Spearman), 
    names_to = "Metric", 
    values_to = "Value"
  ) |>
  ggplot(
    aes(
      x = Value,
      fill = fct_rev(Prediction)
    )
  ) +
  geom_density(
    alpha = 0.5
  ) +
  common_x +
  facet_grid(
    rows = vars(Metric),
    col = vars(Ancestry),
    labeller = labeller(Ancestry = as_labeller(ancestry_labels), 
                        Metric = label_value)
  ) +
  labs(
    x = "Correlation coefficient",
    y = "Density",
    fill = "Prediction"
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Save the plot
ggsave(filename = "Plot_density_dge.pdf", 
       plot = DGE_density_plot, 
       path = path_to_save_location, 
       width = 5, height = 5
       )


# Bar plot ml
ML_bar_plot <- summarized_ml |> 
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
        width = 0.2, position = position_dodge(0.7)
    ) +
    common_y +
    facet_grid(
        rows = vars(Algorithm),
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
      label = ifelse(
        is.na(p_value),  
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p =", format(p_value, digits = 3))  # Otherwise, display the p-value
      )
    ),
    size = 4,    # Adjust text size
    vjust = 2,   # Align text to the top (since we're using Inf for y)
    hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
    inherit.aes = FALSE  # Don't inherit the default aesthetics
  )
# Save 
ggsave(filename = "Plot_bar_ml.pdf", 
       plot = ML_bar_plot, 
       path = path_to_save_location, 
       width = 5, height = 5
       )

# Density plot
ML_density_plot <- metric_ml |>
  pivot_longer(
    cols = c(LogisticRegression, RandomForestClassifier), 
    names_to = "Algorithm", 
    values_to = "Value"
  ) |>
  mutate(
    Value = Value + rnorm(n(), mean = 0, sd = 0.01) # Adding small random noise
  ) |>
  ggplot(
    aes(
      x = Value,
      fill = fct_rev(Prediction)
    )
  ) +
  geom_density(
    alpha = 0.5
  ) +
  common_x +
  facet_grid(
    rows = vars(Algorithm),
    col = vars(Ancestry),
    labeller = labeller(Ancestry = as_labeller(ancestry_labels))
  ) +
  labs(
    x = "ROC AUC",
    y = "Density",
    fill = "Prediction"
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Save the plot
ggsave(filename = "Plot_density_ml.pdf", 
       plot = ML_density_plot, 
       path = path_to_save_location, 
       width = 5, height = 5
       )



# Combine model weights
# Read the files 
combined_regression_weights <- data.frame()
combined_forest_weights <- data.frame()
for (folder in match_vscratch_dir){
    regression_file <- file.path(folder, "Feature_importance_LogisticRegression.csv")
    forest_file <- file.path(folder, "Feature_importance_RandomForestClassifier.csv")

    # Seed information
    seed <- sub(".*_(\\d+)$", "\\1", folder)
    seed <- as.numeric(seed)

    # Regression
    regression_weights <- fread(regression_file) 
    regression_weights$Seed = seed
    combined_regression_weights <- bind_rows(combined_regression_weights, regression_weights) 

    # Forest
    forest_weights <- fread(forest_file) 
    forest_weights$Seed = seed
    combined_forest_weights <- bind_rows(combined_forest_weights, forest_weights) 
}

# Save
fwrite(combined_regression_weights, file.path(path_to_save_location, "Feature_importance_LogisticRegression.csv"))
fwrite(combined_forest_weights, file.path(path_to_save_location, "Feature_importance_RandomForestClassifier.csv"))























