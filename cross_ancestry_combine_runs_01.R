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
p_pearson <- permutation_test(
  contrast_metric_dge, 
  value = "Pearson", 
  group = "Status",
  paired = "Seed",
  tail = "left"
  )
p_spearman <- permutation_test(
  contrast_metric_dge, 
  value = "Spearman", 
  group = "Status",
  paired = "Seed",
  tail = "left"
  )

# ML
p_regression_auc <- filter(contrast_metric_ml, Algorithm == "LogisticRegression") |> 
  permutation_test(
    value = "ROC_AUC", 
    group = "Status",
    paired = "Seed",
    tail = "left"
    )
p_regression_acc <- filter(contrast_metric_ml, Algorithm == "LogisticRegression") |> 
  permutation_test(
    value = "Accuracy", 
    group = "Status",
    paired = "Seed",
    tail = "left"
    )

p_forest_auc <- filter(contrast_metric_ml, Algorithm == "RandomForestClassifier") |>
  permutation_test(
    value = "ROC_AUC", 
    group = "Status",
    paired = "Seed",
    tail = "left"
    )
p_forest_acc <- filter(contrast_metric_ml, Algorithm == "RandomForestClassifier") |>
  permutation_test(
    value = "Accuracy", 
    group = "Status",
    paired = "Seed",
    tail = "left"
    )

# Aggregate data
# DGE - contrast_metric_dge
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
fwrite(aggregated_contrast_metric_dge, 
  file.path(path_to_save_location, "Aggregeated_contrast_metric_dge.csv"))

# DGE - raw_metric_dge
aggregated_raw_metric_dge <- raw_metric_dge |>
  filter(
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

# ML - contrast_metric_ml
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

# ---------------------------------------------------------------------------------------------------------------------
# Baseline
# DGE
baseline_metric_dge <- fload_data(match_vscratch_dir, file = "Baseline_metric_dge.csv")
baseline_metric_per_gene_dge <- fload_data(match_vscratch_dir, file = "Baseline_metric_per_gene.csv")

# Statistics
# DGE
p_cond1_rmse <- filter(baseline_metric_dge, Condition == setup$classification$comparison[1]) |>
  permutation_test(
    value = "RMSE", 
    group = "Status",
    paired = "Seed",
    tail = "right"
    )
p_cond2_rmse <- filter(baseline_metric_dge, Condition == setup$classification$comparison[2]) |>
  permutation_test(
    value = "RMSE", 
    group = "Status",
    paired = "Seed",
    tail = "right"
    )

p_cond1_r2 <- filter(baseline_metric_dge, Condition == setup$classification$comparison[1]) |>
  permutation_test(
    value = "R2", 
    group = "Status",
    paired = "Seed",
    tail = "left"
    )
p_cond2_r2 <- filter(baseline_metric_dge, Condition == setup$classification$comparison[2]) |>
  permutation_test(
    value = "R2", 
    group = "Status",
    paired = "Seed",
    tail = "left"
  )

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
fwrite(aggregated_baseline_metric_dge, 
  file.path(path_to_save_location, "Aggregated_baseline_metric_dge.csv"))

# ---------------------------------------------------------------------------------------------------------------------
# Interactions and baseline 
modified_pattern <- gsub("cross_ancestry", "interactions", match_pattern)
vscratch_dir_in_combined_runs <- file.path("data", "combined_runs")
# Extracting all folders in the 'vscratch_dir_in' that match 'modified_pattern'
modified_all_vscratch_dir_in <- list.dirs(vscratch_dir_in_combined_runs, full.names = TRUE, recursive = FALSE)
modified_match_vscratch_dir <- grep(modified_pattern, modified_all_vscratch_dir_in, value = TRUE)

# Load data
interactions <- fload_data(modified_match_vscratch_dir, file = "Interactions.csv")
baseline <- fload_data(modified_match_vscratch_dir, file = "Baseline.csv")

# Filter sig. features
logFC_threshold <- 1
sig_interactions <- interactions |>
  filter(adj.P.Val < 0.05 & abs(logFC) > logFC_threshold) |>
  pull(Feature)

sig_baseline <- baseline |>
  filter(adj.P.Val < 0.05 & abs(logFC) > logFC_threshold) |>
  pull(Feature)

# Top interactions
top_interactions <- interactions |>
  filter(
    adj.P.Val < 0.05 & 
    abs(logFC) > logFC_threshold
  ) |>
  slice_max(abs(logFC), n = 5) |>
  pull(Feature)


# ---------------------------------------------------------------------------------------------------------------------
# Contrast
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

# Define limits for y-axis
max_y_value <- 1.0 # not considering error bars
max_sd <- 0.1 # estimated max sd
max_multiplicator <- 1.2 + max_sd
# y-limit
max_y_limit <- max_y_value * max_multiplicator
# y-scale
common_y <- scale_y_continuous(
  limits = c(0.0, max_y_limit),  
  breaks = c(0.0, 0.5, 1.0)  
)

# Define common x-axis
common_x <- scale_x_continuous(
  limits = c(0.5, 1.05), 
  breaks = c(0.5, 0.75, 1))

# Plots
# ---- Contrast correlation logFC (DGE) ----
y_max <- max(
  aggregated_contrast_metric_dge$mean_value + aggregated_contrast_metric_dge$sd_value, 
  na.rm = TRUE
  )

# Plot
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
  geom_text(
    aes(
      label = round(mean_value, 3),  
      y = mean_value * 0.5  
    ),
    size = 2,  
    color = "white",  
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - sd_value, 
      ymax = mean_value + sd_value
      ), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  common_y +
  # coord_cartesian(
  #   ylim = c(0.75, 1)
  # ) +
  facet_grid(
    rows = vars(Metric),
    col = vars(Ancestry),
    # labeller = labeller(Ancestry = as_labeller(ancestry_labels))
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  ) +
  # P-value annotation
  # Horizontal line
  annotate(
    "segment",
    x = 1, xend = 2,  
    y = y_max * 1.1, yend = y_max * 1.1,  
    linewidth = 0.3
  ) +
  # Ticks
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = y_max * 1.2,
      label = ifelse(
        is.na(p_value),
        "p[paired~perm.~test] == NA",
        paste0(
          "p[paired~perm.~test] == ", 
          ifelse(p_value < 0.01, 
                format(p_value, digits = 3, scientific = TRUE), 
                format(p_value, digits = 3)),
          ifelse(p_value < 0.001, " **",  # Two asterisks for p < 0.001
                ifelse(p_value < 0.05, " *", ""))  # One asterisk for p < 0.05
        )
      )
    ),
    parse = TRUE,  
    size = 2
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

# ---- Contrast distribution avg. logFC (DGE) ----
# Train - Inference
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
thr = 1
upper_threshold <- predicted_values + thr
lower_threshold <- predicted_values - thr

# Calculate the correlation
cor_pearson <- cor(lm_data$Train, lm_data$Inference, method = "pearson")
cor_p_pearson <- pvalue(independence_test(Inference ~ Train, data = lm_data))

cor_spearman <- cor(lm_data$Train, lm_data$Inference, method = "spearman")
cor_p_spearman <- pvalue(independence_test(rank(Inference) ~ rank(Train), data = lm_data))

# Label
pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

# Add threshold
lm_data <- lm_data |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > thr,
  ) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",                       
      TRUE ~ "Rest"   
      )
  )

# Train and inf ancestry
x_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Train"])
y_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Inference"])

# Compared conditions
condition_1 <- unique(raw_metric_dge$coef)[1]
condition_2 <- unique(raw_metric_dge$coef)[2]

# Legend 
# Legend x-position
max_name_length <- max(nchar(unique(lm_data$coloring)))
legend_x_pos <- min(0.01, 0.01 + 0.005 * max_name_length)  
# Legend y-position
max_items <- length(unique(lm_data$coloring))
legend_y_pos <- max(0.6, 1.0 - 0.01 * max_items)

# Plot
contrast_distribution_logFC <- lm_data |>
  ggplot(
    aes(
      x = Train,
      y = Inference
    )
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = coloring),
    size = point_size,
    show.legend = FALSE
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = "dummy"),
    size = NA,   
    show.legend = FALSE 
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Baseline"),
    aes(color = coloring),
    size = point_size
  ) +
  geom_point(
    data = lm_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
    aes(color = coloring),
    size = point_size
  ) +
  geom_smooth(
    method = "lm",
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +
  geom_line(
    aes(x = Train, y = upper_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) + 
  geom_line(
    aes(x = Train, y = lower_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey",
      "dummy" = "black"
    ),
    breaks = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      "Rest",
      "dummy" 
    ),
    labels = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      pearson_label,
      spearman_label
    )
  ) +
  guides(
    color = guide_legend(
      title = "Ancestry"
      )
  ) +
  labs(
    x = paste(x_ancestry, condition_1, "vs", condition_2, "(avg. logFC)"),
    y = paste(y_ancestry, condition_1, "vs", condition_2, "(avg. logFC)")
  ) +
  geom_text_repel(
    data = lm_data |> filter(Feature %in% top_interactions),
    aes(label = Feature),
    size = 1.5,
    color = "black",  
    segment.color = "black",  
    min.segment.length = 0,
    max.iter = 2000,
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.justification = c(0, 1),
    legend.position = c(legend_x_pos, legend_y_pos), 
    legend.spacing.x = unit(0, "cm"), 
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent"),  
    legend.key.spacing = unit(0, "cm")
  ) 

# Save
ggsave(filename = "Contrast_distribution_avg_logFC_train_inference.pdf", 
       plot = contrast_distribution_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# Train - Subset
lm_data <- aggregated_raw_metric_dge |>
  select(-c(Ancestry, sd_value, se_value)) |>
  pivot_wider(names_from = Status, values_from = mean_value)

# Model
lm_model <- lm(Test ~ Train, data = lm_data)

# Coefficients
predicted_values <- fitted(lm_model)
residuals <- residuals(lm_model)
abs_residuals <- abs(residuals)

# Create threshold line
thr = 1
upper_threshold <- predicted_values + thr
lower_threshold <- predicted_values - thr

# Calculate the correlation
cor_pearson <- cor(lm_data$Train, lm_data$Test, method = "pearson")
cor_p_pearson <- pvalue(independence_test(Test ~ Train, data = lm_data))

cor_spearman <- cor(lm_data$Train, lm_data$Test, method = "spearman")
cor_p_spearman <- pvalue(independence_test(rank(Test) ~ rank(Train), data = lm_data))

# Label
pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

# Add coloring
lm_data <- lm_data |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > thr,
  ) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",                       
      TRUE ~ "Rest"   
      )
  )

# x and y ancestry
x_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Train"])
y_ancestry <- unique(aggregated_raw_metric_dge$Ancestry[aggregated_raw_metric_dge$Status == "Test"])

# Compared conditions
condition_1 <- unique(raw_metric_dge$coef)[1]
condition_2 <- unique(raw_metric_dge$coef)[2]

# Legend 
# Legend x-position
max_name_length <- max(nchar(unique(lm_data$coloring)))
legend_x_pos <- min(0.01, 0.01 + 0.005 * max_name_length)  
# Legend y-position
max_items <- length(unique(lm_data$coloring))
legend_y_pos <- max(0.6, 1.0 - 0.01 * max_items)

# Plot
contrast_distribution_logFC <- lm_data |>
  ggplot(
    aes(
      x = Train,
      y = Test
    )
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = coloring),
    size = point_size,
    show.legend = FALSE
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = "dummy"),
    size = NA,   
    show.legend = FALSE 
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Baseline"),
    aes(color = coloring),
    size = point_size
  ) +
  geom_point(
    data = lm_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
    aes(color = coloring),
    size = point_size
  ) +
  geom_smooth(
    method = "lm",
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +
  geom_line(
    aes(x = Train, y = upper_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) + 
  geom_line(
    aes(x = Train, y = lower_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey",
      "dummy" = "black"
    ),
    breaks = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      "Rest",
      "dummy" 
    ),
    labels = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      pearson_label,
      spearman_label
    )
  ) +
  guides(
    color = guide_legend(
      title = "Subset"
      )
  ) +
  labs(
    x = paste(x_ancestry, condition_1, "vs", condition_2, "(avg. logFC)"),
    y = paste(y_ancestry, condition_1, "vs", condition_2, "(avg. logFC)")
  ) +
  geom_text_repel(
    data = lm_data |> filter(Feature %in% top_interactions),
    aes(label = Feature),
    size = 1.5,
    color = "black",  
    segment.color = "black",  
    min.segment.length = 0,
    max.iter = 2000,
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.justification = c(0, 1),
    legend.position = c(legend_x_pos, legend_y_pos), 
    legend.spacing.x = unit(0, "cm"), 
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent"),  
    legend.key.spacing = unit(0, "cm")
  ) 

# Save
ggsave(filename = "Contrast_distribution_avg_logFC_train_test.pdf", 
       plot = contrast_distribution_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast distribution logFC (DGE) ----
# Train - Inference
lm_data <- raw_metric_dge |>
  filter(coef == unique(raw_metric_dge$coef)[1]) |>
  select(-c(Ancestry)) |>
  pivot_wider(names_from = Status, values_from = logFC) 

# Model
lm_model <- lm(Inference ~ Train, data = lm_data)

# Coefficients
predicted_values <- fitted(lm_model)
residuals <- residuals(lm_model)
abs_residuals <- abs(residuals)

# Create threshold line
thr = 1
upper_threshold <- predicted_values + thr
lower_threshold <- predicted_values - thr

# Calculate the correlation
cor_pearson <- cor(lm_data$Train, lm_data$Inference, method = "pearson")
cor_p_pearson <- pvalue(independence_test(Test ~ Inference, data = lm_data))

cor_spearman <- cor(lm_data$Train, lm_data$Test, method = "spearman")
cor_p_spearman <- pvalue(independence_test(rank(Test) ~ rank(Inference), data = lm_data))

# Label
pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

# Add coloring
lm_data <- lm_data |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > thr,
  ) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",                       
      TRUE ~ "Rest"   
      )
  ) |>
  mutate(
    Seed = factor(Seed)
  )

# x and y ancestry
x_ancestry <- unique(raw_metric_dge$Ancestry[raw_metric_dge$Status == "Train"])
y_ancestry <- unique(raw_metric_dge$Ancestry[raw_metric_dge$Status == "Inference"])

# Compared conditions
condition_1 <- unique(raw_metric_dge$coef)[1]
condition_2 <- unique(raw_metric_dge$coef)[2]

# Plot
contrast_distribution_logFC <- lm_data |>
  ggplot(
    aes(
      x = Train,
      y = Inference
    )
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(
      color = coloring,
      shape = Seed,
      ),
    size = point_size,
    show.legend = FALSE
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = "dummy"),
    size = NA,   
    show.legend = FALSE 
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Baseline"),
    aes(
      color = coloring,
      shape = Seed
      ),
    size = point_size
  ) +
  geom_point(
    data = lm_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
    aes(
      color = coloring,
      shape = Seed
      ),
    size = point_size
  ) +
  geom_smooth(
    method = "lm",
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +
  geom_line(
    aes(x = Train, y = upper_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) + 
  geom_line(
    aes(x = Train, y = lower_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey",
      "dummy" = "black"
    ),
    breaks = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      "Rest",
      "dummy" 
    ),
    labels = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      pearson_label,
      spearman_label
    )
  ) +
  guides(
    color = guide_legend(
      title = "Ancestry",
      title.theme = element_text(face = "bold")
      ),  
    shape = guide_legend(
      title = "Seed",
      ncol = 5,
      title.theme = element_text(face = "plain")
      )
  ) +
  labs(
    x = paste(x_ancestry, condition_1, "vs", condition_2, "(logFC)"),
    y = paste(y_ancestry, condition_1, "vs", condition_2, "(logFC)")
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    # legend.title = element_blank(),
    legend.justification = c(0, 1),
    legend.position = c(legend_x_pos, legend_y_pos), 
    legend.spacing.x = unit(0, "cm"), 
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent"),  
    legend.key.spacing = unit(0, "cm")
  ) 

# Save
ggsave(filename = "Contrast_distribution_logFC_train_inference.pdf", 
       plot = contrast_distribution_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# Train - Subset
lm_data <- raw_metric_dge |>
  filter(coef == unique(raw_metric_dge$coef)[1]) |>
  select(-c(Ancestry)) |>
  pivot_wider(names_from = Status, values_from = logFC) 

# Model
lm_model <- lm(Test ~ Train, data = lm_data)

# Coefficients
predicted_values <- fitted(lm_model)
residuals <- residuals(lm_model)
abs_residuals <- abs(residuals)

# Create threshold line
thr = 1
upper_threshold <- predicted_values + thr
lower_threshold <- predicted_values - thr

# Calculate the correlation
cor_pearson <- cor(lm_data$Train, lm_data$Test, method = "pearson")
cor_p_pearson <- pvalue(independence_test(Test ~ Train, data = lm_data))

cor_spearman <- cor(lm_data$Train, lm_data$Test, method = "spearman")
cor_p_spearman <- pvalue(independence_test(rank(Test) ~ rank(Train), data = lm_data))

# Label
pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

# Add coloring
lm_data <- lm_data |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > thr,
  ) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",                       
      TRUE ~ "Rest"   
      )
  ) |>
  mutate(
    Seed = factor(Seed)
  )

# x and y ancestry
x_ancestry <- unique(raw_metric_dge$Ancestry[raw_metric_dge$Status == "Train"])
y_ancestry <- unique(raw_metric_dge$Ancestry[raw_metric_dge$Status == "Test"])

# Compared conditions
condition_1 <- unique(raw_metric_dge$coef)[1]
condition_2 <- unique(raw_metric_dge$coef)[2]

# Plot
contrast_distribution_logFC <- lm_data |>
  ggplot(
    aes(
      x = Train,
      y = Test
    )
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(
      color = coloring,
      shape = Seed,
      ),
    size = point_size,
    show.legend = FALSE
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Rest"),
    aes(color = "dummy"),
    size = NA,   
    show.legend = FALSE 
  ) +
  geom_point(
    data = lm_data |> filter(coloring == "Baseline"),
    aes(
      color = coloring,
      shape = Seed
      ),
    size = point_size
  ) +
  geom_point(
    data = lm_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
    aes(
      color = coloring,
      shape = Seed
      ),
    size = point_size
  ) +
  geom_smooth(
    method = "lm",
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +
  geom_line(
    aes(x = Train, y = upper_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) + 
  geom_line(
    aes(x = Train, y = lower_threshold), 
    color = "blue", 
    linetype = "dashed",
    linewidth = point_size,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey",
      "dummy" = "black"
    ),
    breaks = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      "Rest",
      "dummy" 
    ),
    labels = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      pearson_label,
      spearman_label
    )
  ) +
  guides(
    color = guide_legend(
      title = "Subset",
      title.theme = element_text(face = "bold")
      ),  
    shape = guide_legend(
      title = "Seed",
      ncol = 5,
      title.theme = element_text(face = "plain")
      )
  ) +
  labs(
    x = paste(x_ancestry, condition_1, "vs", condition_2, "(logFC)"),
    y = paste(y_ancestry, condition_1, "vs", condition_2, "(logFC)")
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    # legend.title = element_blank(),
    legend.justification = c(0, 1),
    legend.position = c(legend_x_pos, legend_y_pos), 
    legend.spacing.x = unit(0, "cm"), 
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent"),  
    legend.key.spacing = unit(0, "cm")
  ) 

# Save
ggsave(filename = "Contrast_distribution_logFC_train_test.pdf", 
       plot = contrast_distribution_logFC, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast density R-values (DGE) ----
contrast_density_R <- contrast_metric_dge |>
  pivot_longer(
    cols = c(
      Pearson, 
      Spearman
      ), 
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
    color = NA,
    alpha = 0.5
  ) +
  common_x +
  facet_grid(
    rows = vars(Metric),
    col = vars(Ancestry),
    # labeller = labeller(Ancestry = as_labeller(ancestry_labels))
  ) +
  labs(
    x = "Correlation coefficient",
    y = "Density",
    fill = "Prediction"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.15, 0.9),  
    legend.direction = "vertical"
  ) 

# Save
ggsave(filename = "Contrast_density_R_values.pdf", 
       plot = contrast_density_R, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast prediction phenotype (ML) ----
y_max <- max(
  aggregated_contrast_metric_ml$mean_value + aggregated_contrast_metric_ml$sd_value, 
  na.rm = TRUE
  )

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
  annotate(
    "segment",
    x = 1, xend = 2,  
    y = y_max * 1.1, yend = y_max * 1.1,  
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = y_max * 1.2,
      label = ifelse(
        is.na(p_value),
        "p[paired~perm.~test] == NA",
        paste0(
          "p[paired~perm.~test] == ", 
          ifelse(p_value < 0.01, 
                format(p_value, digits = 3, scientific = TRUE), 
                format(p_value, digits = 3)),
          ifelse(p_value < 0.001, " **",  # Two asterisks for p < 0.001
                ifelse(p_value < 0.05, " *", ""))  # One asterisk for p < 0.05
        )
      )
    ),
    parse = TRUE,  
    size = 2
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
threshold <- 0.5
# Conditions
class_0 <- colnames(raw_metric_ml)[1]
class_1 <- colnames(raw_metric_ml)[2] # Target

contrast_distribution_probabilities <- raw_metric_ml |>
  mutate(
    Classification = ifelse(.data[[class_1]] > threshold, class_1, class_0),
    True_labels = as.factor(y)
  ) |>
  select(-all_of(class_0)) |>
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
    linewidth = point_size
  ) + 
  facet_grid(
    cols = vars(fct_rev(Prediction)),
    rows = vars(Algorithm),
    labeller = labeller(
      Algorithm = as_labeller(algorithm_labels)
    )
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
  theme_small_legend() +
  theme(
    legend.position = c(0.9, 0.1)
  )

# Save
ggsave(filename = "Contrast_distribution_probabilities.pdf", 
       plot = contrast_distribution_probabilities, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Contrast density prediction (ML) -----
contrast_density_prediction <- contrast_metric_ml |>
  pivot_longer(
    cols = c(ROC_AUC, Accuracy), 
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
    color = NA,
    alpha = 0.5
  ) +
  common_x +
  facet_grid(
    rows = vars(Metric),
    col = vars(Ancestry, Algorithm),
    labeller = labeller(
      # Ancestry = as_labeller(ancestry_labels),
      Algorithm = as_labeller(algorithm_labels),
      Metric = as_labeller(metric_labels)
      )
  ) +
  labs(
    x = "X",
    y = "Density",
    fill = "Prediction"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.15, 0.9),  
    legend.direction = "vertical"
  ) 

# Save
ggsave(filename = "Contrast_density_prediction_values.pdf", 
       plot = contrast_density_prediction, 
       path = path_to_save_location, 
       width = 3, height = 3
       )


# ---------------------------------------------------------------------------------------------------------------------
# Baseline
# Visualize
condition_labels <- aggregated_baseline_metric_dge |>
  distinct(Condition, n_condition) |>
  mutate(label = paste0(Condition, " (n = ", n_condition, ")")) |>
  pull(label, Condition) 

# Define limits for y-axis
max_y_value <- 1.5 # not considering error bars
max_sd <- 0.1 # estimated max sd
max_multiplicator <- 1.2 + max_sd
# y-limit
max_y_limit <- max_y_value * max_multiplicator
# y-scale
common_y <- scale_y_continuous(
  limits = c(0.0, max_y_limit),  
  breaks = c(0.0, 0.5, 1.0, 1.5)  
)

# Plots
# ---- Prediction baseline (DGE) ----
y_max <- max(
  aggregated_baseline_metric_dge$mean_value + aggregated_baseline_metric_dge$sd_value, 
  na.rm = TRUE
  )

# Plot
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
    # P-value annotation
  # Horizontal line
  annotate(
    "segment",
    x = 1, xend = 2,  
    y = y_max * 1.1, yend = y_max * 1.1,  
    linewidth = 0.3
  ) +
  # Ticks
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_max * 1.08, yend = y_max * 1.105,
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = y_max * 1.2,
      label = ifelse(
        is.na(p_value),
        "p[paired~perm.~test] == NA",
        paste0(
          "p[paired~perm.~test] == ", 
          ifelse(p_value < 0.01, 
                format(p_value, digits = 3, scientific = TRUE), 
                format(p_value, digits = 3)),
          ifelse(p_value < 0.001, " **",  # Two asterisks for p < 0.001
                ifelse(p_value < 0.05, " *", ""))  # One asterisk for p < 0.05
        )
      )
    ),
    parse = TRUE,  
    size = 2
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

# ---- Highest error genes (DGE) ----
genes <- baseline_metric_per_gene_dge |>
  filter(Prediction == "Ancestry") |>
  group_by(
    Condition,
    Feature
  ) |>                
  summarize(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Condition, desc(mean_RMSE)) |>
  group_by(Condition) |>    
  slice_max(mean_RMSE, n = 5, with_ties = FALSE) |>
  select(Condition, Feature) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",                       
      TRUE ~ "Rest"   
      )
  )

# Plot
highest_error_genes <- baseline_metric_per_gene_dge |>
  inner_join(genes, by = c("Condition", "Feature")) |>
  ggplot(
    aes(
      x = Feature,
      y = RMSE,
      color = coloring
    )
  ) +
  geom_boxplot(
    outlier.size = 0.1
  ) +
  facet_grid(
    cols = vars(fct_rev(Prediction)),
    rows = vars(Condition),
    scales = "free_x"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3) 
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey"),
    breaks = c("Interaction + Baseline", "Interaction", "Baseline"),
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(filename = "Baseline_highest_error_genes.pdf", 
       plot = highest_error_genes, 
       path = path_to_save_location, 
       width = 3, height = 3
       )