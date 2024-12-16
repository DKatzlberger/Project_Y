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
  # When the script is run from the command line then 'output_directory' is given
  # The pattern to extract all matchin directories is extracted from 'output_directory'
  output_path = setup$output_directory
  # This regular expression removes trailing underscores followed by digits from strings 
  match_pattern <- gsub("_\\d+$", "", output_path)
  # TODO - Check if regular expression also works for other approaches

} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "job_settings.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

  # When run for development no 'output_directory' is specified
  # Hence, the 'output_directory' has to be contructed like in the 'process_yaml_... .py'

  # Construction:
  # Vscratch_dir is the place where the files are stored
  vscratch_dir_in = "data/runs"
  # Tag is used to specify which data it is e.g. TCGA, NIAGADS
  tag <- setup$tag
  # Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
  comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
  # Analysis_name specifies what was analysed
  # E.g. comparsion of ancestries EUR_to_AMR, subsetting_EUR
  # This often is modified depending which analysis
  train_ancestry <- toupper(setup$classification$train_ancestry)
  infer_ancestry <- toupper(setup$classification$infer_ancestry)
  analysis_name  <-  paste0(train_ancestry, "_to_", infer_ancestry)
  # Combination of components to create the 'match_pattern'
  # The 'match_pattern' is used as pattern to extract all folders in the vscratch dir
  match_pattern <- paste0(comparison, "_", analysis_name)
}

# Extracting all folders in the 'vscratch_dir' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vsratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Check if there were matching folders
if (length(match_vsratch_dir) == 0) {
  message("No matching folders found.")
}

# Save the results of the analysis
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- "data/combined_runs"
path_to_save_location <- file.path(vscratch_dir_out, comparison, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}


# Combining metric csv files 
# Each folder is one seed
metric_dge <- data.frame()
metric_ml <- data.frame()
for (folder in matching_folders){
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data for each seed
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data for each seed
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

# metric_dge |>
#   select(!c(V1, V2, Comparison, n_ancestry))

# Add otlier column to check if ROC_AUC values contain outliers
metric_ml <- metric_ml |>
  group_by(Status) |>
  mutate(
    ROC_AUC_outlier = hampel_filter(ROC_AUC),  # Hampel filter for ROC AUC
  )

# Differential gene expression analysis:

# Test difference in distribution of correlated logFC  
# 1. Kolmogorow-Smirnow
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

# 3. t-test
# parametric
pearson_test <- metric_dge |>
  filter(Status == "Test",
         #Seed %in% c(1, 2, 3, 4, 5)
         ) |>
  pull(Pearson)
pearson_inference <- metric_dge |>
  filter(Status == "Inference",
         #Seed %in% c(1, 2, 3, 4, 5)
         ) |>
  pull(Pearson)

t_pearson <- t.test(pearson_test, pearson_inference)

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
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
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
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
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



# Save DGE and ML approach
# 1. Metric
fwrite(metric_dge, file.path(path_to_save_location, "Metric_dge.csv"))
fwrite(metric_ml, file.path(path_to_save_location, "Metric_ml.csv"))
# 2. Summarized metric
fwrite(summarized_dge, file.path(path_to_save_location, "Summarized_metric_dge.csv"))
fwrite(summarized_ml, file.path(path_to_save_location, "Summarized_metric_ml.csv"))
# 3. Plots
# Save the image
ggsave(filename = "Plot_dge.pdf", plot = DGE_plot, 
       path = path_to_save_location, width = 5, height = 5
       )

ggsave(filename = "Plot_ml.pdf", plot = ML_plot, 
       path = path_to_save_location, width = 5, height = 5
       )

# Save the R object
saveRDS(DGE_plot, file = file.path(path_to_save_location, "Plot_dge.rds"))
saveRDS(ML_plot, file = file.path(path_to_save_location, "Plot_ml.rds"))


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
fwrite(combined_weights, file.path(path_to_save_location, "Weights.csv"))

# TODO - Create overall directory in runs to move single run files






















