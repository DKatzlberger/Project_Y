# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Visualization
    library(patchwork)
    library(ggrepel)
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
  yaml_file <- "PanCanAtlas_RSEM_covariates_subtype_age_EUR_to_EAS.yml"
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
  analysis_name  <-  "EUR_subsetting"
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

# Save the results of the analysis: 'EUR_subsetting'
# 'Vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- "data/combined_runs"
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}


# Loading and combining the matched folders
metric_dge <- data.frame()
metric_ml <- data.frame()
for (folder in match_vsratch_dir){
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data for each seed
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data for each seed
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}

# Summarize metric across seeds
# Calculate mean, sd, sem 
# 1. DGE
# 2. ML
summarized_dge <- metric_dge |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = "Value",
        names_to = "Metric"
    ) |> 
    # Dont use 'Seed' in 'group_by'
    group_by(Metric, 
             n_test_ancestry,
    ) |>  
    summarize(
        n_seeds = n(),
        mean_value = mean(Value, na.rm = TRUE),
        sd_value = sd(Value, na.rm = TRUE),
        se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
  ) |>
  mutate(Algorithm = "limma")

summarized_ml <- metric_ml |> 
    pivot_longer(
        cols = c(LogisticRegression, RandomForestClassifier),
        values_to = "Value",
        names_to = "Algorithm"
    ) |> 
    # Dont use 'Seed' in 'group_by'
    group_by(Algorithm, 
             n_test_ancestry,
             Metric
    ) |>  
    summarize(
        n_seeds = n(),
        mean_value = mean(Value, na.rm = TRUE),
        sd_value = sd(Value, na.rm = TRUE),
        se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
  )

# Combine DGE and ML 
# Add 'Metric_Type' column
combined_metric <- bind_rows(summarized_dge, summarized_ml) |>
  mutate(Metric_type = paste0(Metric, " (", Algorithm, ")")) |>
    mutate(
    Metric_type = factor(Metric_type, 
                         levels = c(
                           "ROC_AUC (LogisticRegression)", 
                           "ROC_AUC (RandomForestClassifier)", 
                           "Spearman (limma)",
                           "Pearson (limma)"
                         ))
  )

# Save the combined metric
fwrite(combined_metric, file.path(path_to_save_location, "Combined_metric.csv"))

# Visualize
performance_plot <- combined_metric |>
  ggplot(
    aes(
      x = n_test_ancestry,
      y = mean_value,
      color = Metric_type
    )
  ) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_value - sd_value, 
      ymax = mean_value + sd_value
      ), 
    width = 0.1
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = c(0.0, 0.5, 1.0)
  ) +
  labs(
    #title = "Performance per sample size",
    x = "Test sample size (EUR)",
    y = "Mean value",
    color = "Method"
  )

train_number_plot <- combined_metric |>
  head(1) |>
  ggplot(
    aes(
      x = Ancestry,
      y = n_train_ancestry
    )
  ) +
  geom_col() +
  labs(
    x = "Train ancestry",
    y = "Sample size"
  )

# Combine plots
combined_plot <- train_number_plot + performance_plot + plot_layout(widths = c(1, 4)) +
  plot_annotation(title = "Performance per sample size")
 
# Save
ggsave(filename = "Performance_by_sample_size.pdf", 
       plot = performance_plot, 
       path = path_to_save_location, 
       width = 10, height = 5
       )

