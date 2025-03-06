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
    library(limma)
    library(fgsea)
    # Visualization
    library(patchwork)
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

# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Save the results of the analysis
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}


# Contrast -----------------------------------------------------------------------------------
# Load  data 
metric_dge <- fload_data(
  folders = match_vscratch_dir,
  file = "Metric_contrast.csv"
)

metric_ml <- fload_data(
  folders = match_vscratch_dir,
  file = "Metric_ml.csv"
)

# Gene metric
contrast_leading_edge <- fload_data(
  folders = match_vscratch_dir,
  file = "Metric_contrast_per_gene.csv"
)

# Save 
fwrite(metric_dge, file.path(path_to_save_location, "Metric_dge.csv"))
fwrite(metric_ml, file.path(path_to_save_location, "Metric_ml.csv"))
fwrite(contrast_leading_edge, file.path(path_to_save_location, "Contrast_metric_per_gene.csv"))


# Differential gene expression -------------------------------------------------------------------------
# Permutation testing (Returns NA if no variation across subsets)
# Non-parametric
# Without assumption
p_pearson <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Pearson")
p_spearman <- permutation_test(data = metric_dge, group_col = "Status", value_col = "Spearman")

p_rmse <- permutation_test(data = metric_dge, group_col = "Status", value_col = "RMSE")
p_r2 <- permutation_test(data = metric_dge, group_col = "Status", value_col = "R2")

# Summarize data (taking mean of correlation values)
summarized_dge <- metric_dge |> 
    pivot_longer(
        cols = c(
          Pearson, 
          Spearman,
          RMSE,
          R2
          ),
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
  mutate(
    p_value = case_when(
      Metric == "Pearson" ~ p_pearson,
      Metric == "Spearman" ~ p_spearman,
      Metric == "RMSE" ~ p_rmse,
      Metric == "R2" ~ p_r2
    )
  )

# Machine learning -------------------------------------------------------------------------
# Permutation testing
# Non-parametric
# Without assumption
p_regression_auc <- metric_ml |>
  filter(Algorithm == "LogisticRegression") |>
  permutation_test(group_col = "Status", value_col = "ROC_AUC")

p_regression_acc <- metric_ml |>
  filter(Algorithm == "LogisticRegression") |>
  permutation_test(group_col = "Status", value_col = "Accuracy")


p_forest_auc <- metric_ml |>
  filter(Algorithm == "RandomForestClassifier") |>
  permutation_test(group_col = "Status", value_col = "ROC_AUC")

p_forest_acc <- metric_ml |>
  filter(Algorithm == "RandomForestClassifier") |>
  permutation_test(group_col = "Status", value_col = "Accuracy")

# Prepare data format
summarized_ml <- metric_ml |> 
    pivot_longer(
        cols = c(ROC_AUC, Accuracy),
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
fwrite(summarized_dge, file.path(path_to_save_location, "Summarized_metric_dge.csv"))
fwrite(summarized_ml, file.path(path_to_save_location, "Summarized_metric_ml.csv"))

# Visualization -------------------------------------------------------------------------------------
# Common settings for plots 
# Labeller to annotate number of samples per ancestry
ancestry_labels <- summarized_dge |>
  distinct(Ancestry, n_inf_ancestry) |>
  mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
  select(Ancestry, label) |>
  deframe()

# Axis
common_y <- scale_y_continuous(
  limits = c(0, 1.5), 
  breaks = c(0, 0.5, 1))

common_x <- scale_x_continuous(
  limits = c(0, 1.2), 
  breaks = c(0, 0.5, 1))

# Differential analysis --------------------------------------------------------------------------------------------
# ---- Correlation of logFC ----
# Bar
correlation_bar_plot <- summarized_dge |> 
  filter(Metric %in% c("Pearson", "Spearman")) |>
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
    title = "Correlation of logFC",
    x = "Prediction",
    y = "Correlation coefficient"
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p = ", 
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
ggsave(filename = "Plot_bar_correlation_of_logFC.pdf", 
       plot = correlation_bar_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )
# ---- Density ----
correlation_density_plot <- metric_dge |>
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
    labeller = labeller(Ancestry = as_labeller(ancestry_labels))
  ) +
  labs(
    title = "Distribution of correlation coefficients",
    x = "Correlation coefficient",
    y = "Density",
    fill = "Prediction"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# Save the plot
ggsave(filename = "Plot_density_correlation_of_logFC.pdf", 
       plot = correlation_density_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Predicting logFC (prediction_bar_plot) -----
# Bar 
prediction_bar_plot <- summarized_dge |> 
  filter(Metric %in% c("RMSE", "R2")) |>
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
    title = "Prediction of logFC",
    x = "Prediction",
    y = "Y"
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p = ", 
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
ggsave(filename = "Plot_bar_prediction_of_logFC.pdf", 
       plot = prediction_bar_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Leading edge (contrast_leading_edge_plot) ----
# Features are select based on highest error in ancestry
leading_edge_genes <- contrast_leading_edge |>
  filter(Prediction == "Ancestry") |>
  group_by(
    Comparison,
    Feature
  ) |>                
  summarize(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Comparison, desc(mean_RMSE)) |>
  group_by(Comparison) |>    
  slice_max(mean_RMSE, n = 20) |>
  pull(Feature) |>
  unique()

# Plot
contrast_leading_edge_plot <- contrast_leading_edge |>
  filter(
    Feature %in% leading_edge_genes
  ) |>
  ggplot(
    aes(
      x = Feature,
      y = RMSE
    )
  ) +
  geom_boxplot(
    outlier.size = 0.1
  ) +
  facet_grid(
    cols = vars(Ancestry, fct_rev(Prediction)),
    labeller = labeller(
      Ancestry = as_labeller(ancestry_labels)
      ),
    scales = "free_x"
  ) +
  labs(
    title = "Features with biggest error in ancestry"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(
     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 10))
  )

# Save 
ggsave(filename = "Plot_boxplot_leading_edge_contrast.pdf", 
       plot = contrast_leading_edge_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Patchwork (patchwork_prediction_of_contrast) ----
patchwork_prediction_of_contrast <- prediction_bar_plot + 
  contrast_leading_edge_plot

# Save
ggsave(filename = "Patchwork_prediction_of_contrast.pdf", 
       plot = patchwork_prediction_of_contrast, 
       path = path_to_save_location, 
       width = 6, height = 3
       )


# Comparison of gene errors
# Limma
# Step 1: Reshape Data - Keep all seeds as observations
# error_long <- contrast_leading_edge |>
#   filter(Comparison == contrast_leading_edge$Comparison[1]) |>
#   select(Feature, Status, RMSE, Seed)

# # Step 2: Convert to matrix (genes as rows, samples as columns)
# error_matrix <- error_long |>
#   pivot_wider(names_from = c(Seed, Status), values_from = RMSE) |>
#   column_to_rownames("Feature") |>
#   as.matrix()

# # Step 3: Create meta
# observations  <- colnames(error_matrix)
# status <- ifelse(grepl("Test", colnames(error_matrix)), "Test", "Inference")
# meta <- data.frame(
#   Observation = observations, 
#   Status = status) |>
#   mutate(Status = factor(Status, levels = c("Test", "Inference")))
# # Check if matrix and meta align
# meta <- meta[match(colnames(rmse_matrix), meta$Observation), ] 

# # Step 4: Create design matrix for limma
# design <- model.matrix(~Status, data = meta)

# # Step 5: Fit the model
# fit <- lmFit(error_matrix, design)
# fit <- eBayes(fit)
# contrast_predicted_gene_errors <- extract_results(fit)

# Permutation test
# Step 1: Reshape Data - Keep all seeds as observations
error_long <- contrast_leading_edge |>
  filter(Comparison == contrast_leading_edge$Comparison[1]) |>
  select(Feature, Status, RMSE, Seed) |>
  mutate(Status = factor(Status, levels = c("Test", "Inference"))) |>
  as.data.table()

# Step 2: Perform permutation test
process_batch_parallel <- function(batch_data) {
  # Calculate the p-values from the independence test
  p_values <- batch_data[, .(p_value = as.numeric(pvalue(independence_test(
    RMSE ~ Status, 
    data = .SD
  ))[[1]])), by = Feature]
  
  # Calculate the mean RMSE for each Feature and Status group
  mean_rmse <- batch_data[, .(mean_rmse = mean(RMSE, na.rm = TRUE)), by = .(Feature, Status)]
  
  # Reshape the data to calculate the log2FC between 'Test' and 'Inference'
  mean_rmse_wide <- dcast(mean_rmse, Feature ~ Status, value.var = "mean_rmse")
  
  # Calculate log2 fold change (log2FC)
  mean_rmse_wide[, logFC := log2(`Inference` / `Test`)]  # log2 of FC
  
  # Merge p-values with log2FC
  result <- merge(p_values, mean_rmse_wide[, .(Feature, logFC)], by = "Feature")
  
  return(result)
}

# Set up the number of cores (adjust to your machine's capabilities)
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export required objects to each worker
clusterExport(cl, list("error_long", "process_batch_parallel", "independence_test"))

# Load required libraries on each worker
clusterEvalQ(cl, {
  library(coin)
  library(data.table)
})

# Split data into batches (100 features per batch)
batch_size <- 100
feature_list <- unique(error_long$Feature)
batches <- split(feature_list, ceiling(seq_along(feature_list) / batch_size))

# Use parLapply to process each batch in parallel
perm_results <- do.call(rbind, parLapply(cl, batches, function(batch) {
  # Subset the data for the current batch
  batch_data <- error_long[Feature %in% batch]
  # Perform the permutation test for the batch
  process_batch_parallel(batch_data)
}))

# Stop the cluster when done
stopCluster(cl)

# Step 3: Adjust p_values for multiple testing
perm_results[, p_adjusted := p.adjust(p_value, method = "fdr")]

# Add ancestry
perm_results[, Ancestry := unique(contrast_leading_edge$Ancestry)]

# Save
fwrite(perm_results, file.path(path_to_save_location, "Contrast_metric_errorFC.csv"))

# Functional analysis
# Fgsea
database <- "data/downloads/geneset-libraries/MSigDB_Hallmark_2020.txt"
database <- read_enrichR_database(database)
enrichment <- perform_gsea(perm_results, database, rank_by = "logFC")
enrichment[, Ancestry := infer_ancestry]

# Save
fwrite(enrichment, file.path(path_to_save_location, "Contrast_enrichment.csv"))

# Visualize -----------------------------------------------------------------------------------------
# ---- Volcano_plot (contrast_gene_errors_volcano_plot) ----
logFC_threshold <- 1
point_size  <- 0.5
contrast_gene_errors_volcano_plot <- perm_results |>
   ggplot(
    aes(
      x = logFC,
      y = -log10(p_adjusted),
      color = (p_adjusted < 0.05 & abs(logFC) > logFC_threshold)
    )
  ) +
  geom_point(size = point_size) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "grey"),
  ) +
  labs(
    title = "Volcano plot of gene errors",
    y = "-log10(adj.P.Val)"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.position = "none"
  )

# Save
ggsave(filename = "Contrast_gene_errors_volcano_plot.pdf", 
       plot = contrast_gene_errors_volcano_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Leading egdge (contrast_leading_edge_gene_errors_plot) ----
contrast_leading_edge_gene_errors_plot <- perm_results |>
  filter(p_adjusted < 0.05 & abs(logFC) > logFC_threshold) |>
  slice_max(logFC, n = 30) |>
  ggplot(
    aes(
      x = Feature,
      y = logFC
    )
  ) +
  geom_col(
    width = 0.7
  ) +
  labs(
    title = "Features with highest logFC in ancestry",
    x = "Feature",
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 10))
  )

# Save
ggsave(filename = "Contrast_leading_edge_gene_errors_bar_plot.pdf", 
       plot = contrast_leading_edge_gene_errors_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Patchwork ----
patchwork_contrast_logFC_gene_error <- contrast_gene_errors_volcano_plot +
  contrast_leading_edge_gene_errors_plot

# Save
ggsave(filename = "Patchwork_contrast_logFC_gene_error.pdf", 
       plot = patchwork_contrast_logFC_gene_error, 
       path = path_to_save_location, 
       width = 6, height = 3
       )


# Machine learning ----------------------------------------------------------------------------------
# Custom labeller 
algorithm_labels <- c(
  "LogisticRegression" = "Logistic Regression",
  "RandomForestClassifier" = "Random Forest")

metric_labels <- c(
  "ROC_AUC" = "ROC AUC",
  "Accuracy" = "Accuracy"
)

# --- Predicting phenotype ----
prediction_phenotype_bar_plot <- summarized_ml |> 
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
    col = vars(Ancestry, Algorithm),
    labeller = labeller(Ancestry = as_labeller(ancestry_labels),
                        Algorithm = as_labeller(algorithm_labels),
                        Metric = as_labeller(metric_labels)
                      )
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p = ", 
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
    title = "Prediction of phenotypes",
    x = "Prediction",
    y = "Y"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip()

# Save 
ggsave(filename = "Plot_bar_prediction_of_phenotypes.pdf", 
       plot = prediction_phenotype_bar_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Density ----
prediction_phenotypes_density_plot <- metric_ml |>
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
      Ancestry = as_labeller(ancestry_labels),
      Algorithm = as_labeller(algorithm_labels),
      Metric = as_labeller(metric_labels)
      )
  ) +
  labs(
    title = "Distribution of prediction scores",
    x = "ROC AUC",
    y = "Density",
    fill = "Prediction"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
  ) 

# Save 
ggsave(filename = "Plot_density_prediction_of_phenotypes.pdf", 
       plot = prediction_phenotypes_density_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )


# Baseline ------------------------------------------------------------------------------------------
# Load data
baseline_metric <- fload_data(
  folders = match_vscratch_dir,
  file = "Metric_baseline.csv"
)

# Statisitical test (perm_test)
p_cond1_rmse <- baseline_metric |>
  filter(Condition == setup$classification$comparison[1]) |>
  permutation_test(group_col = "Status", value_col = "RMSE")

p_cond2_rmse <- baseline_metric |>
  filter(Condition == setup$classification$comparison[2]) |>
  permutation_test(group_col = "Status", value_col = "RMSE")


p_cond1_r2 <- baseline_metric |>
  filter(Condition == setup$classification$comparison[1]) |>
  permutation_test(group_col = "Status", value_col = "R2")

p_cond2_r2 <- baseline_metric |>
  filter(Condition == setup$classification$comparison[2]) |>
  permutation_test(group_col = "Status", value_col = "R2")

# Summarize
summarized_baseline_metric <- baseline_metric |>
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
fwrite(summarized_baseline_metric, file.path(path_to_save_location, "Baseline_summarized_metric_dge.csv"))


baseline_leading_edge <- data.frame()
for (folder in match_vscratch_dir){
    file <- file.path(folder, "Metric_baseline_per_gene.csv")

    # Load and append DGE data for each seed
    data <- fread(file) 
    baseline_leading_edge <- bind_rows(baseline_leading_edge, data) 
}
# Save
fwrite(baseline_leading_edge, file.path(path_to_save_location, "Baseline_metric_per_gene.csv"))

# Visualize ------------------------------------------------------------------------------------------
condition_labels <- summarized_baseline_metric |>
  distinct(Condition, n_condition) |>
  mutate(label = paste0(Condition, " (n = ", n_condition, ")")) |>
  pull(label, Condition) 

# ---- Prediction of baseline (prediction_of_baseline_plot) ----
prediction_of_baseline_plot <- summarized_baseline_metric |>
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
    col = vars(Ancestry, Condition),
    labeller = labeller(
      Ancestry = as_labeller(ancestry_labels),
      Condition = as_labeller(condition_labels)
      )
  ) +
  geom_text(
    aes(
      x = 0.5,  # Align text to the left side of the first bar
      y = Inf,  # Position the text at the top of the plot
      label = ifelse(
        is.na(p_value),
        paste("Perm. test,", "p = NA"),  
        paste("Perm. test,", "p = ", 
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
    title = "Prediction of baseline",
    x = "Prediction",
    y = "Y"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip () 

# Save
ggsave(filename = "Plot_bar_prediction_of_baseline.pdf", 
       plot = prediction_of_baseline_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Leading edge (baseline_leading_edge_plot) ----
leading_edge_genes <- baseline_leading_edge |>
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
  slice_max(mean_RMSE, n = 20, with_ties = FALSE) |>
  select(Condition, Feature)

# Plot
baseline_leading_edge_plot <- baseline_leading_edge |>
  inner_join(leading_edge_genes, by = c("Condition", "Feature")) |>
  # filter(
  #   # Prediction == "Ancestry",
  #   Feature %in% leading_edge_genes
  # ) |>
  ggplot(
    aes(
      x = Feature,
      y = RMSE
    )
  ) +
  geom_boxplot(
    outlier.size = 0.1
  ) +
  facet_grid(
    cols = vars(Ancestry, Condition),
    rows = vars(fct_rev(Prediction)),
    labeller = labeller(
      Ancestry = as_labeller(ancestry_labels),
      Condition = as_labeller(condition_labels)
      ),
    scales = "free_x"
  ) +
  labs(
    title = "Features with biggest error in ancestry"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(
     axis.text.x = element_text(angle = 90, hjust = 1)
  )

# Save
ggsave(filename = "Plot_boxplot_leading_edge_baseline.pdf", 
       plot = baseline_leading_edge_plot, 
       path = path_to_save_location, 
       width = 3, height = 3
       )

# ---- Patchwork (patchwork_prediction_of_baseline) ----
patchwork_prediction_of_baseline <- prediction_of_baseline_plot + 
  baseline_leading_edge_plot

# Save
ggsave(filename = "Patchwork_prediction_of_baseline.pdf", 
       plot = patchwork_prediction_of_baseline, 
       path = path_to_save_location, 
       width = 6, height = 3
       )



# Machine learning interpretations -----------------------------------------------------------------------------------------------
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

# Similarity of the EUR-subset ---------------------------------------------------------------------------------------------------
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





















