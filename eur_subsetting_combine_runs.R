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
source("figure_themes.R")

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
  yaml_file <- "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_EAS.yml"
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
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Save the results of the analysis: 'EUR_subsetting'
vscratch_dir_out  <- "data/combined_runs"
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Contrast
# DGE
contrast_metric_dge <- fload_data(match_vscratch_dir, file = "Contrast_metric_dge.csv")
raw_metric_dge <- fload_data(match_vscratch_dir, file = "LogFCs.csv")

# ML
contrast_metric_ml <- fload_data(match_vscratch_dir, file = "Contrast_metric_ml.csv")
raw_metric_ml <- fload_data(match_vscratch_dir, file = "Probabilities.csv")

# Aggregate data
group_seeds <- c(2, 5, 8, "all")
# DGE - contrast_metric_dge
aggregated_dfs <- list()
for (seed_count in group_seeds) {
  # If seed_count is "all", use all the seeds in the data
  if (seed_count == "all") {
    filtered_data <- contrast_metric_dge
  } else {
    # Randomly sample the seeds based on the specified seed_count
    set.seed(42)  # For reproducibility
    sampled_seeds <- sample(unique(contrast_metric_dge$Seed), seed_count)
    
    # Filter the data for the randomly sampled seeds
    filtered_data <- contrast_metric_dge |>
      filter(Seed %in% sampled_seeds)
  }

  aggregated_df <- filtered_data |> 
    pivot_longer(
      cols = c(Pearson, Spearman),
      values_to = "Value",
      names_to = "Metric"
    ) |> 
    group_by(
      Metric, 
      n_test_ancestry,
    ) |>  
    summarize(
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(
      Algorithm = "Limma",
      Method = "Statistics",
      Aggregation = "R-values",
      Metric_type = paste0(Metric, " (", Aggregation, ")")
    )

  # Store the aggregated dataframe in the list
  aggregated_dfs[[paste("aggregated_seed_", seed_count, sep = "")]] <- aggregated_df
}
aggregated_contrast_metric_dge <- bind_rows(aggregated_dfs)

# DGE - raw_metric_dge
aggregated_raw_dfs <- list()
for (seed_count in group_seeds) {
  
  # If seed_count is "all", use all the seeds in the data
  if (seed_count == "all") {
    filtered_data <- raw_metric_dge
  } else {
    # Randomly sample the seeds based on the specified seed_count
    set.seed(123)  # For reproducibility
    sampled_seeds <- sample(unique(raw_metric_dge$Seed), seed_count)
    
    # Filter the data for the randomly sampled seeds
    filtered_data <- raw_metric_dge |>
      filter(Seed %in% sampled_seeds)
  }
  
  # Perform aggregation for the filtered data
  aggregated_raw_df <- filtered_data |>
    filter(coef == unique(filtered_data$coef)[1]) |>
    group_by(coef, Feature, n_test_ancestry, Ancestry) |>
    summarise(
      n_seeds = n(),
      mean_value = mean(logFC, na.rm = TRUE),
      sd_value = sd(logFC, na.rm = TRUE),
      se_value = sd(logFC, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Correlation
  wide <- aggregated_raw_df |>
    select(n_test_ancestry, Feature, mean_value) |>
    pivot_wider(names_from = n_test_ancestry, values_from = mean_value)
  
  cor_pearson <- cor(wide[,-1], method = "pearson")
  cor_spearman <- cor(wide[,-1], method = "spearman")
  
  # Correlation dataframes
  cor_pearson_df <- tibble(
      n_test_ancestry = rownames(cor_pearson), 
      mean_value = cor_pearson[, "constant"],
      sd_value = NA,
      n_seeds = unique(aggregated_raw_df$n_seeds)
    ) |>
    filter(n_test_ancestry != "constant") |>
    mutate(n_test_ancestry = as.numeric(n_test_ancestry)) |>
    mutate(
      Algorithm = "Limma",
      Method = "Statistics",
      Aggregation = "logFCs",
      Metric = "Pearson",
      Metric_type = paste0(Metric, " (", Aggregation, ")")
    )
  
  cor_spearman_df <- tibble(
      n_test_ancestry = rownames(cor_spearman), 
      mean_value = cor_spearman[, "constant"], 
      sd_value = NA,
      n_seeds = unique(aggregated_raw_df$n_seeds)
    )  |>
    filter(n_test_ancestry != "constant") |>
    mutate(n_test_ancestry = as.numeric(n_test_ancestry)) |>
    mutate(
      Algorithm = "Limma",
      Method = "Statistics",
      Aggregation = "logFCs",
      Metric = "Spearman",
      Metric_type = paste0(Metric, " (", Aggregation, ")")
    )
  
  # Combine across metrics
  combined_contrast_raw_metric <- bind_rows(
    cor_pearson_df, 
    cor_spearman_df
  )
  
  # Store the aggregated dataframe in the list
  aggregated_raw_dfs[[paste("aggregated_seed_", seed_count, sep = "")]] <- combined_contrast_raw_metric
}
aggregated_raw_metric_dge <- bind_rows(aggregated_raw_dfs)

# ML - contrast_metric_ml
aggregated_dfs <- list()
for (seed_count in group_seeds) {
  
  # If "all" is specified, use all available seeds
  if (seed_count == "all") {
    filtered_data <- contrast_metric_ml
  } else {
    # Set seed for reproducibility
    set.seed(123)  
    # Randomly sample the required number of seeds
    sampled_seeds <- sample(unique(contrast_metric_ml$Seed), seed_count)
    
    # Filter data based on sampled seeds
    filtered_data <- contrast_metric_ml |> 
      filter(Seed %in% sampled_seeds)
  }
  
  # Perform aggregation
  aggregated_ml_df <- filtered_data |> 
    pivot_longer(
      cols = c(ROC_AUC, Accuracy),
      values_to = "Value",
      names_to = "Metric"
    ) |> 
    group_by(
      Algorithm, 
      n_test_ancestry,
      Metric
    ) |>  
    summarize(
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(
      Method = "Machine learning",
      Metric_type = paste0(Metric, " (", Algorithm, ")")
    ) |>
    mutate(
      Metric_type = recode(
        Metric_type,
        "Accuracy (LogisticRegression)" = "Accuracy (Logistic Regression)",
        "ROC_AUC (LogisticRegression)" = "ROC AUC (Logistic Regression)",
        "Accuracy (RandomForestClassifier)" = "Accuracy (Random Forest)",
        "ROC_AUC (RandomForestClassifier)" = "ROC AUC (Random Forest)"
      )
    )
  
  # Store in list
  aggregated_dfs[[paste("aggregated_seed_", seed_count, sep = "")]] <- aggregated_ml_df
}
aggregated_contrast_metric_ml <- bind_rows(aggregated_dfs)

# Visualize
point_size = 0.5

algorithm_labels <- c(
  "LogisticRegression" = "Logistic Regression",
  "RandomForestClassifier" = "Random Forest")

metric_labels <- c(
  "ROC_AUC" = "ROC AUC",
  "Accuracy" = "Accuracy"
)


# ---- Perfromance statisitc (DGE) ----
# Combine across approaches
combined_contrast_metric_dge <- bind_rows(
  aggregated_contrast_metric_dge,
  aggregated_raw_metric_dge
  ) 

performance_plot <- combined_contrast_metric_dge |>
  ggplot(
    aes(
      x = n_test_ancestry,
      y = mean_value,
      color = Metric,
    )
  ) +
  geom_point(
    size = point_size
  ) +
  geom_line(
    linewidth = (point_size / 2)
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - sd_value, 
      ymax = mean_value + sd_value
      ), 
    linewidth = (point_size / 2)
  ) +
  facet_grid(
    cols = vars(Aggregation),
    rows = vars(n_seeds)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 3) 
  ) +
  labs(
    x = "n (Testset)",
    y = "Y",
    color = "Method"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    panel.grid.major = element_line(color = "lightgrey", size = (point_size / 3)),
    legend.title = element_blank(),
    legend.position = c(1, 0.01),  
    legend.direction = "vertical",
    legend.justification = c(1, 0)
  ) 

# Save
ggsave(filename = "Performance_by_sample_size_statistic.pdf", 
       plot = performance_plot, 
       path = path_to_save_location, 
       width = 5, height = 3
       )
# ---- Perfromance machine learning (ML) ----
performance_plot_1 <- aggregated_contrast_metric_ml |>
  ggplot(
    aes(
      x = n_test_ancestry,
      y = mean_value,
      color = Metric
    )
  ) +
  geom_point(
    size = point_size
  ) +
  geom_line(
    linewidth = (point_size / 2)
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - sd_value, 
      ymax = mean_value + sd_value
      ), 
    width = 0.1
  ) +
  facet_grid(
    cols = vars(Algorithm),
    rows = vars(n_seeds),
    labeller = labeller(
      Algorithm = as_labeller(algorithm_labels)
      )
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = c(0.0, 0.5, 1.0)
  ) +
  scale_color_discrete(
    labels = metric_labels
  ) +
  labs(
    x = "n (Testset)",
    y = "Y",
    color = "Method"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    panel.grid.major = element_line(color = "lightgrey", size = (point_size / 3)),
    legend.title = element_blank(),
    legend.position = c(1, 0.01),  
    legend.direction = "vertical",
    legend.justification = c(1, 0)
  ) 

# Save
ggsave(filename = "Performance_by_sample_size_machine_learning.pdf", 
       plot = performance_plot_1, 
       path = path_to_save_location, 
       width = 4, height = 3
       )

# ---- Distriubution probabilities (ML) ----
threshold <- 0.5
# Conditions
class_0 <- colnames(raw_metric_ml)[1]
class_1 <- colnames(raw_metric_ml)[2] 

distribution_probabilities <- raw_metric_ml |>
  select(-all_of(class_0)) |>
  mutate(
    True_labels = as.factor(y)
  ) |>
  pivot_longer(
    cols = 1,
    names_to = "Condition",
    values_to = "Probabilities"
  ) |>
  ggplot(
    aes(
      x = Condition,
      y = Probabilities,
      fill = True_labels
    )
  ) +
  geom_violin(
    color = NA
  ) +
  geom_hline(
    yintercept = threshold, 
    color = "blue", 
    linetype = "dashed", 
    linewidth = point_size / 2
  ) + 
  facet_grid(
    cols = vars(factor(n_test_ancestry)),
    rows = vars(Algorithm),
    labeller = labeller(
      Algorithm = as_labeller(algorithm_labels)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0.0, 0.5, 1.0)
  ) +
  scale_fill_discrete(
    name = "True labels",
    labels = c(class_0, class_1)
  ) +
  labs(
    x = "Target condition"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() 

# Save
ggsave(filename = "Distribution_probabilities.pdf", 
       plot = distribution_probabilities, 
       path = path_to_save_location, 
       width = 5, height = 3
       )






