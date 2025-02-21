# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Parallelization
    library(parallel)
    # Statistics 
    library(coin)
    library(fgsea)
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
  yaml_file <- "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_ADMIX.yml"
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

# Save location
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# ---------------------------------------------------------------------------------------------------------------------
# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Contrast
# DGE
contrast_metric_dge <- fload_data(match_vscratch_dir, file = "Contrast_metric_dge.csv")
raw_metric_dge <- fload_data(match_vscratch_dir, file = "LogFCs.csv")

# ML
contrast_metric_ml <- fload_data(match_vscratch_dir, file = "Contrast_metric_ml.csv")
raw_metric_ml <- fload_data(match_vscratch_dir, file = "Probabilities.csv")

# Statisitcs
# DGE
p_pearson <- permutation_test(contrast_metric_dge, "Pearson", "Status")
p_spearman <- permutation_test(contrast_metric_dge, "Spearman", "Status")

# ML
p_regression_auc <- filter(contrast_metric_ml, Algorithm == "LogisticRegression") |> 
  permutation_test("ROC_AUC", "Status")
p_regression_acc <- filter(contrast_metric_ml, Algorithm == "LogisticRegression") |> 
  permutation_test("Accuracy", "Status")

p_forest_auc <- filter(contrast_metric_ml, Algorithm == "RandomForestClassifier") |>
  permutation_test("ROC_AUC", "Status")
p_forest_acc <- filter(contrast_metric_ml, Algorithm == "RandomForestClassifier") |>
  permutation_test("Accuracy", "Status")

# Aggregate data
# DGE - metric_contrast_dge
aggregated_contrast_metric_dge <- contrast_metric_dge |>
  pivot_longer(
    cols = c(
      Pearson, 
      Spearman,
      # RMSE,
      # R2
      ),
    values_to = "Value",
    names_to = "Metric"
  ) |> 
  group_by(
    Metric, 
    Ancestry, 
    Status, 
    Prediction, 
    n_inf_ancestry
  ) |>  
  summarize(
    n_seeds = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    p_value = case_when(
      Metric == "Pearson" ~ p_pearson,
      Metric == "Spearman" ~ p_spearman,
      # Metric == "RMSE" ~ p_rmse,
      # Metric == "R2" ~ p_r2
    )
  )

# Save
fwrite(aggregeated_contrast_metric_dge, 
  file.path(path_to_save_location, "Aggregeated_contrast_metric_dge.csv"))

# DGE - metric_raw_dge
aggregated_raw_metric_dge <- raw_metric_dge |>
  filter(
    Status %in% c("Train", "Inference"),
    coef == unique(raw_metric_dge$coef)[1]
  ) |>
  group_by(
    coef,
    Feature,
    Status,
    Ancestry,
    n_train_ancestry,
    n_inf_ancestry,
  ) |>
  summarise(
    n_seeds = n(),
    mean_value = mean(logFC, na.rm = TRUE),
    sd_value = sd(logFC, na.rm = TRUE),
    se_value = sd(logFC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Save
fwrite(aggregated_raw_metric_dge, 
  file.path(path_to_save_location, "Aggregeated_raw_metric_dge.csv"))

# ML
aggregated_contrast_metric_ml <- contrast_metric_ml |>
 pivot_longer(
    cols = c(
      ROC_AUC, 
      Accuracy
      ),
    values_to = "Value",
    names_to = "Metric"
  ) |> 
  group_by(
    Metric,
    Algorithm,
    Ancestry, 
    Status, 
    Prediction, 
    n_inf_ancestry
  ) |>  
  summarize(
    n_seeds = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
      p_value = case_when(
        Metric == "ROC_AUC" & Algorithm == "LogisticRegression" ~ p_regression_auc,
        Metric == "Accuracy" & Algorithm == "LogisticRegression" ~ p_regression_acc,
        Metric == "ROC_AUC" & Algorithm == "RandomForestClassifier" ~ p_forest_auc,
        Metric == "Accuracy" & Algorithm == "RandomForestClassifier" ~ p_forest_acc
      )
  )

# Save
fwrite(aggregated_contrast_metric_ml, 
  file.path(path_to_save_location, "Aggregated_contrast_metric_ml.csv"))


# Visualize
point_size = 0.5

ancestry_labels <- aggregated_contrast_metric_ml |>
  distinct(Ancestry, n_inf_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()

algorithm_labels <- c(
  "LogisticRegression" = "Logistic Regression",
  "RandomForestClassifier" = "Random Forest")

metric_labels <- c(
  "ROC_AUC" = "ROC AUC",
  "Accuracy" = "Accuracy"
)

common_y <- scale_y_continuous(
  limits = c(0.0, 1.5), 
  breaks = c(0.0, 0.5, 1))

# ---- Contrast correlation logFC (DGE) ----
contrast_correlation_logFC <- aggregated_contrast_metric_dge |>
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
    labeller = labeller(Ancestry = as_labeller(ancestry_labels))
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("p = NA"),  
        paste("p = ", 
            ifelse(p_value < 0.01, 
                   format(p_value, digits = 3, scientific = TRUE), 
                   format(p_value, digits = 3)))
                   )
      ),
      size = 1.5,    # Adjust text size
      vjust = 1.5,   # Align text to the top
      hjust = 0,     # Align text to the left
      inherit.aes = FALSE  # Don't inherit the default aesthetics
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip()

# Save
ggsave(filename = "Contrast_correlation_logFC.pdf", 
       plot = contrast_correlation_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast distribution logFC (DGE) ----
lm_data <- aggregated_raw_metric_dge |>
  select(-c(Ancestry, sd_value, se_value)) |>
  pivot_wider(names_from = Status, values_from = mean_value)

# Model
lm_model <- lm(Inference ~ Train, data = lm_data)

# Coefficients
predicted_values <- fitted(lm_model)
residuals <- residuals(lm_model)
abs_residuals <- abs(residuals)

# Create threshold line
thr = 2
upper_threshold <- predicted_values + thr
lower_threshold <- predicted_values - thr

# Add threshold
lm_data <- lm_data |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > thr 
  ) 

x_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Train"])
y_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Inference"])

contrast_distribution_logFC <- lm_data |>
  ggplot(
    aes(
      x = Train,
      y = Inference
    )
  ) +
  geom_point(
    aes(color = above_threshold),
    size = point_size
  ) +
  geom_smooth(
    method = "lm",
    linewidth = point_size
  ) +
  geom_line(
    aes(x = Train, y = upper_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.7
  ) + 
  geom_line(
    aes(x = Train, y = lower_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("grey", "red")
  ) +
  geom_text_repel(
    data = lm_data %>% filter(above_threshold), 
    aes(label = Feature), 
    size = 1.5,  
    max.overlaps = 10  
  ) +
  labs(
    x = paste("LogFC", x_ancestry),
    y = paste("LogFC", y_ancestry)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(legend.position = "none")

# Save
ggsave(filename = "Contrast_distribution_logFC.pdf", 
       plot = contrast_distribution_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast prediction phenotype (ML) ----
contrast_prediction_phenotype <- aggregated_contrast_metric_ml |>
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
    cols = vars(Algorithm),
    rows = vars(Metric),
    labeller = labeller(
      Algorithm = as_labeller(algorithm_labels),
      Metric = as_labeller(metric_labels)
      )
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("p = NA"),  
        paste("p = ", 
            ifelse(p_value < 0.01, 
                   format(p_value, digits = 3, scientific = TRUE), 
                   format(p_value, digits = 3)))
                   )
      ),
      size = 1.5,    # Adjust text size
      vjust = 1.5,   # Align text to the top
      hjust = 0,     # Align text to the left
      inherit.aes = FALSE  # Don't inherit the default aesthetics
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip()

# Save
ggsave(filename = "Contrast_prediction_phenotype.pdf", 
       plot = contrast_prediction_phenotype, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast distribution probabilities (ML) ----
contrast_distribution_probabilities <- raw_metric_ml |>
  pivot_longer(
    cols = 1:2,
    names_to = "Condition",
    values_to = "Probabilities"
  ) |>
  ggplot(
    aes(
      x = Condition,
      y = Probabilities
    )
  ) +
  geom_violin() +
  geom_hline(
    yintercept = 0.5, 
    color = "blue", 
    linetype = "dashed", 
    linewidth = point_size
  ) + 
  facet_grid(
    cols = vars(fct_rev(Prediction)),
    rows = vars(Algorithm),
    labeller = labeller(
      Algorithm = as_labeller(algorithm_labels)
    )
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_background() +
  theme_white_strip()

# Save
ggsave(filename = "Contrast_distribution_probabilities.pdf", 
       plot = contrast_distribution_probabilities, 
       path = path_to_save_location, 
       width = 3, height = 3
       )


# Baseline
# DGE
baseline_metric_dge <- fload_data(match_vscratch_dir, file = "Baseline_metric_dge.csv")

# Statistics
# DGE
p_cond1_rmse <- filter(baseline_metric_dge, Condition == setup$classification$comparison[1]) |>
  permutation_test("RMSE", "Status")

p_cond2_rmse <- filter(baseline_metric_dge, Condition == setup$classification$comparison[2]) |>
  permutation_test("RMSE", "Status")


p_cond1_r2 <- filter(baseline_metric_dge, Condition == setup$classification$comparison[1]) |>
  permutation_test("R2", "Status")

p_cond2_r2 <- filter(baseline_metric_dge, Condition == setup$classification$comparison[2]) |>
  permutation_test("R2", "Status")

# Aggregate data
# DGE - baseline_metric_dge
aggregated_baseline_metric_dge <- baseline_metric_dge |>
  pivot_longer(
    cols = c(
      RMSE,
      R2
      ),
    values_to = "Value",
    names_to = "Metric"
  ) |> 
  group_by(
    Condition,
    Ancestry,
    Status, 
    Prediction, 
    Metric, 
    n_inf_ancestry,
    n_condition
  ) |>  
  summarize(
    n_seeds = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    p_value = case_when(
      Metric == "RMSE" & Condition == setup$classification$comparison[1] ~ p_cond1_rmse,
      Metric == "RMSE" & Condition == setup$classification$comparison[2] ~ p_cond2_rmse,
      Metric == "R2" & Condition == setup$classification$comparison[1] ~ p_cond1_r2,
      Metric == "R2" & Condition == setup$classification$comparison[2] ~ p_cond2_r2
    )
  )

# Save
fwrite(aggregated_baseline_metric_dge, file.path(path_to_save_location, "Aggregated_baseline_metric_dge.csv"))

# Visualize
condition_labels <- aggregated_baseline_metric_dge |>
  distinct(Condition, n_condition) |>
  mutate(label = paste0(Condition, " (n = ", n_condition, ")")) |>
  pull(label, Condition) 

# ---- Prediction baseline (DGE) ----
baseline_prediction_condition <- aggregated_baseline_metric_dge |>
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
    col = vars(Condition),
    labeller = labeller(
      Condition = as_labeller(condition_labels)
      )
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("p = NA"),  
        paste("p = ", 
            ifelse(p_value < 0.01, 
                   format(p_value, digits = 3, scientific = TRUE), 
                   format(p_value, digits = 3)))
                   )
      ),
      size = 1.5,    # Adjust text size
      vjust = 1.5,   # Align text to the top
      hjust = 0,     # Align text to the left
      inherit.aes = FALSE  # Don't inherit the default aesthetics
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() 

# Save
ggsave(filename = "Baseline_prediction_condition.pdf", 
       plot = baseline_prediction_condition, 
       path = path_to_save_location, 
       width = 3, height = 3
       )