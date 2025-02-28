# Remove start up messages
suppressPackageStartupMessages(
  {
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # Statistics library
    library(edgeR)
    # Visualization
    library(patchwork)
  }
)
# Custom functions
source("r_utils.R")

# Here starts the script
print("Evoking R...")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  yaml_file <- args[1]
  # Check if it's a valid YAML file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
} else {
  print("Running interactive mode for development.")
  # Dev settings if no command-line argument provided
  yaml_file <- "dev_settings.yml"
  setup <- yaml.load_file(yaml_file)
}

# Transform settings into R useable form
comparison <- setup$classification$comparison
output_column <- setup$classification$output_column
ancestry_column <- setup$classification$ancestry_column
train_ancestry <- setup$classification$train_ancestry
infer_ancestry <- setup$classification$infer_ancestry

# Set seed (because why not)
set.seed(setup$seed)

# Load data
adata <- read_h5ad(setup$data_path)

# Load the observation used in ml
train_idx <- yaml.load_file(file.path(setup$output_directory, "Obs_train.yml"))
test_idx <- yaml.load_file(file.path(setup$output_directory, "Obs_test.yml"))
inf_idx <- yaml.load_file(file.path(setup$output_directory, "Obs_inf.yml"))

# Get number of samples
train_n = length(train_idx)
test_n = length(test_idx)
inf_n = length(inf_idx)

# Subset the data by the indexes create in python
train_data <- adata[train_idx, ]
test_data <- adata[test_idx, ]
inf_data <- adata[inf_idx, ]

# Check if saved 'inf_idx' are from the correct ancestry
inf_adata = adata[adata$obs[ancestry_column] == infer_ancestry]
inf_adata = inf_adata[inf_adata$obs[[output_column]] %in% comparison]
inf_adata_obs = inf_adata$obs_names
stopifnot(length(inf_adata_obs) == length(inf_idx))

# Number of samples per condition
condition_count <- group_by(inf_adata$obs, !!sym(output_column)) |> count()

# Covariate:
# 1. Extract the covariate, if it exists and has a meaningful value
covariate_list <- get_covariate(setup)
# 2. Check if there are replicates and at least two levels for the covariate 
if (!is.null(covariate_list)) {
  # Covariate processing
  print("Checking covariates.")
  covariate_types <- classify_covariates(train_data$obs, covariate_list)

  # Check conditions for "Discrete" covariates
  check_covariate_conditions(train_data$obs, covariates = covariate_types$discrete)
  check_covariate_conditions(test_data$obs, covariates = covariate_types$discrete)
  check_covariate_conditions(inf_data$obs, covariates = covariate_types$discrete)

  # Visualization when covariates exist
  if (setup$visual_val) {
    # Visualize all discrete columns
    to_visualize_columns <- c(setup$classification$output_column, covariate_types$discrete)
    # Datasets to visualize
    datasets <- list(
      Train = train_data$obs,
      Test = test_data$obs,
      `Inf` = inf_data$obs
    )
    # Create and save the plot
    patchwork_plot <- create_stratification_plot(
      datasets = datasets,
      to_visualize_columns = to_visualize_columns
    )
    ggsave(
      filename = "Validation_stratification.pdf",
      plot = patchwork_plot,
      path = setup$output_directory,
      width = 5, height = 5
    )
  }
} else {
  # Visualization when no covariates but visual_val is TRUE
  if (setup$visual_val) {
    # Visualize only the output column
    to_visualize_columns <- setup$classification$output_column
    # Datasets to visualize
    datasets <- list(
      Train = train_data$obs,
      Test = test_data$obs,
      `Inf` = inf_data$obs
    )
    # Create and save the plot
    patchwork_plot <- create_stratification_plot(
      datasets = datasets,
      to_visualize_columns = to_visualize_columns
    )
    ggsave(
      filename = "Validation_stratification_no_covariates.pdf",
      plot = patchwork_plot,
      path = setup$output_directory,
      width = 5, height = 5
    )
  }
}

# Design matrix for GLM
output_column <- setup$classification$output_column
train_design <- create_design(output_column, train_data$obs, covariate = covariate_list)                                                       
test_design <- create_design(output_column, test_data$obs, covariate = covariate_list)
inf_design <- create_design(output_column, inf_data$obs, covariate = covariate_list)

# Filtering genes
if (setup$filter_features & setup$data_type == "expression"){
  print("Filter features")
  # Based on train data samples
  dge <- DGEList(counts = t(train_data$X))
  dge <- calcNormFactors(dge)
  cpm_values <- cpm(dge)
  # Filter genes based on CPM threshold 
  # (e.g., keep genes with CPM > 1 in at least 50% of the samples)
  threshold <- 1 
  filter_genes <- rowSums(cpm_values > threshold) >= (0.5 * ncol(cpm_values))
  # Subset features
  train_filtered <- train_data[, filter_genes]
  test_filtered <- test_data[, filter_genes]
  inf_filtered <- inf_data[, filter_genes]
} else{
  train_filtered <- train_data
  test_filtered <- test_data
  inf_filtered <- inf_data
}
# Save feature (for use in ML)
output_directory <- setup$output_directory
write_yaml(train_filtered$var_names, file.path(output_directory, "Features.yml"))


# Limma workflow
print("Start differential gene expression analysis.")
# Select normalization method
data_type <- setup$data_type
dge_normalization <- setup$dge_normalization
normalization_method <- normalization_methods[[data_type]][[dge_normalization]]
# Normalization
train_norm <- normalization_method(train_filtered$X, train_design)
test_norm <- normalization_method(test_filtered$X, test_design)
inf_norm <- normalization_method(inf_filtered$X, inf_design)

# Fit the model (Means model)
train_limma_fit <- lmFit(train_norm, design = train_design)
test_limma_fit <- lmFit(test_norm, design = test_design)
inf_limma_fit <- lmFit(inf_norm, design = inf_design)

# Ebayes
train_limma_fit <- eBayes(train_limma_fit)
test_limma_fit <- eBayes(test_limma_fit)
inf_limma_fit <- eBayes(inf_limma_fit)

# Results means model
train_mean_res <- extract_results(train_limma_fit)
test_mean_res <- extract_results(test_limma_fit)
inf_mean_res <- extract_results(inf_limma_fit)

# Save results
output_directory <- setup$output_directory
fwrite(train_mean_res, file.path(output_directory, "Limma_means_train.csv"))
fwrite(test_mean_res, file.path(output_directory, "Limma_means_test.csv"))
fwrite(inf_mean_res, file.path(output_directory, "Limma_means_inf.csv"))

# Create contrast matrix 
print("Fit contrast.")
# Used for hypothesis testing between groups
comparison <- setup$classification$comparison
contrast_matrix_train <- create_contrast(colnames(train_design), conditions = comparison)
contrast_matrix_test <- create_contrast(colnames(test_design), conditions = comparison)
contrast_matrix_inf <- create_contrast(colnames(inf_design), conditions = comparison)

# Fit contrast
train_limma_fit_contrast <- contrasts.fit(train_limma_fit, contrast_matrix_train)
test_limma_fit_contrast <- contrasts.fit(test_limma_fit, contrast_matrix_test)
inf_limma_fit_contrast <- contrasts.fit(inf_limma_fit, contrast_matrix_inf)

# Ebayes
train_limma_fit_contrast <- eBayes(train_limma_fit_contrast)
test_limma_fit_contrast <- eBayes(test_limma_fit_contrast)
inf_limma_fit_contrast <- eBayes(inf_limma_fit_contrast)

# Results contrast
train_contrast_res <- extract_results(train_limma_fit_contrast)
test_contrast_res <- extract_results(test_limma_fit_contrast)
inf_contrast_res <- extract_results(inf_limma_fit_contrast)

# Save results
output_directory <- setup$output_directory
fwrite(train_contrast_res, file.path(output_directory, "Limma_contrast_train.csv"))
fwrite(test_contrast_res, file.path(output_directory, "Limma_contrast_test.csv"))
fwrite(inf_contrast_res,file.path(output_directory, "Limma_contrast_inf.csv"))

# Filter coef that are in comparison (only interested in those)
train_contrast_res <- filter(train_contrast_res, coef %in% setup$classification$comparison)
test_contrast_res <- filter(test_contrast_res, coef %in% setup$classification$comparison)
inf_contrast_res <- filter(inf_contrast_res, coef %in% setup$classification$comparison)


# Visualization (Validation of models) -----------------------------------------------------------------------------------------
# Check if logFC of models align with signal in the data
if (setup$visual_val) {
  # Get genes for validation
  top <- train_mean_res |>
    slice_max(logFC, n = 5) |>
    pull(Feature) |>
    unique()

  low <- train_mean_res |>
    slice_min(logFC, n = 5) |>
    pull(Feature) |>
    unique()

  # Combine genes of interest
  goi <- c(top, low)

  # Validation of the training
  t <- visual_validation(
    meta = train_data$obs,
    signal = train_norm$E,
    mean_stats = train_mean_res,
    contrast_stats = train_contrast_res,
    goi = goi,
    data_output_column = setup$classification$output_column
  )

  i <- visual_validation(
    meta = inf_data$obs,
    signal = inf_norm$E,
    mean_stats = inf_mean_res,
    contrast_stats = inf_contrast_res,
    goi = goi,
    data_output_column = setup$classification$output_column
  )

  # Save the pictures
  ggsave(file.path(setup$output_directory, "Validation_stats_train.pdf"),
         plot = t, height = 6, width = 12)
  ggsave(file.path(setup$output_directory, "Validation_stats_inf.pdf"),
         plot = i, height = 6, width = 12)
}


# Raw logFC
train_log_FC <- train_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Train",
    Ancestry = toupper(setup$classification$train_ancestry),
  )

test_log_FC <- test_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Test",
    Ancestry = toupper(setup$classification$train_ancestry),
  )

inf_log_FC <- inf_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Inference",
    Ancestry = toupper(setup$classification$infer_ancestry)
  )

# Combine
raw_logFC <- bind_rows(train_log_FC, test_log_FC, inf_log_FC) |>
  mutate(
    Seed = setup$seed,
    n_train_ancestry = train_n,
    n_inf_ancestry = inf_n,
  )

# Save
fwrite(raw_logFC, file.path(setup$output_directory, "LogFCs.csv"))

# Calculate metric
print("Correlation of logFCs.")
# Select and rename logFC columns from each data frame for clarity
train_log_FC <- train_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_train = logFC)

test_log_FC <- test_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_test = logFC)

inf_log_FC <- inf_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_inf = logFC)

# Merge all three data frames (aligned by coef and feature)
merged_log_FC <- train_log_FC |>
  inner_join(test_log_FC, by = c("coef", "Feature")) |>
  inner_join(inf_log_FC, by = c("coef", "Feature"))

# Correlation summarized across genes
pearson <- head(compute_correlation(data = merged_log_FC, method = "pearson"), 2)
spearman <- head(compute_correlation(data = merged_log_FC, method = "spearman"), 2)

# Make metric dataframe
metric_df <- inner_join(pearson, spearman, by = c("V1", "V2")) |>
  mutate(
    Seed = setup$seed,
    Status = ifelse(V1 == "logFC_test", "Test", "Inference"),
    Ancestry = toupper(setup$classification$infer_ancestry),
    Prediction = ifelse(V1 == "logFC_test", "Subset", "Ancestry"),
    n_train_ancestry = train_n,
    n_inf_ancestry = inf_n,
  )

# Save
fwrite(metric_df, file.path(setup$output_directory, "Contrast_metric_dge.csv"))



# # Predicting logFC ---------------------------------------------------------------------------------------------
# print("Predicting logFC.")

# # Predicted contrast (train set)
# predicted_contrast_train <- train_contrast_res |> 
#   select(coef, Feature, logFC) |>
#   rename(Comparison = coef, Predicted_contrast = logFC) |>
#   arrange(Comparison, Feature)

# # Contrast (test set, inf set)
# # Observed response
# meta <- as_tibble(test_filtered$obs, rownames = "patient_id")  |> select(patient_id, setup$classification$output_column)
# observed_response_test <- calculate_observed_mean_response(as_tibble(t(test_norm$E), rownames = "patient_id"), 
#                                                            meta = meta,
#                                                            output_column = setup$classification$output_column)

# meta <- as_tibble(inf_filtered$obs, rownames = "patient_id") |> select(patient_id, setup$classification$output_column)
# observed_response_inf <- calculate_observed_mean_response(as_tibble(t(inf_norm$E), rownames = "patient_id"), 
#                                                           meta = meta,
#                                                           output_column = setup$classification$output_column)

# # Observed contrast
# observed_contrast_test <- calculate_observed_contrast(observed_response_test, 
#                                                       output_column = setup$classification$output_column)

# observed_contrast_inf <- calculate_observed_contrast(observed_response_inf, 
#                                                      output_column = setup$classification$output_column) 

# # ---- Per gene metric ----
# # Test
# test_gene_error <- observed_contrast_test |>
#   left_join(predicted_contrast_train, by = c("Comparison", "Feature")) |>
#   group_by(Comparison, Feature) |>
#   summarise(
#     RMSE = sqrt(mean((Observed_contrast - Predicted_contrast)^2)),  
#     # MAE  = mean(abs(Observed_contrast - Predicted_contrast)),
#     # R2 = 1 - sum((Observed_contrast - Predicted_contrast)^2) / sum((Observed_contrast - mean(Observed_contrast))^2),
#     .groups = "drop"
#   ) |>
#   mutate(
#     Status = "Test",
#     Prediction = "Subset"
#   )

# # Inf
# inf_gene_error <- observed_contrast_inf |>
#   left_join(predicted_contrast_train, by = c("Comparison", "Feature")) |>
#   group_by(Comparison, Feature) |>
#   summarise(
#     RMSE = sqrt(mean((Observed_contrast - Predicted_contrast)^2)),  
#     # MAE  = mean(abs(Observed_contrast - Predicted_contrast)),
#     # R2 = 1 - sum((Observed_contrast - Predicted_contrast)^2) / sum((Observed_contrast - mean(Observed_contrast))^2),
#     .groups = "drop"
#   ) |>
#   mutate(
#     Status = "Inference",
#     Prediction = "Ancestry"
#   )

# # Combine
# contrast_per_gene_metric <- bind_rows(test_gene_error, inf_gene_error) |>
#   mutate(
#     Seed = setup$seed,
#     Ancestry = toupper(setup$classification$infer_ancestry),
#     n_train_ancestry = train_n,
#     n_inf_ancestry = inf_n,
#   )

# # Save
# fwrite(contrast_per_gene_metric, file.path(setup$output_directory, "Metric_contrast_per_gene.csv"))

# # ---- Summarized metric ----
# test_prediction <- observed_contrast_test |>
#   left_join(predicted_contrast_train, by = c("Comparison", "Feature")) |>
#   group_by(Comparison) |>
#   summarise(
#     # Calculate the MAE
#     MAE = mean(abs(Observed_contrast - Predicted_contrast)),
#     # Calculate the RMSE
#     RMSE = sqrt(mean((Observed_contrast - Predicted_contrast)^2)),
#     # Calculate the R-squared
#     SS_residual = sum((Observed_contrast - Predicted_contrast)^2),  
#     SS_total = sum((Observed_contrast - mean(Observed_contrast))^2), 
#     R2 = 1 - (SS_residual / SS_total),
#     .groups = "drop"
#   ) |>
#   mutate(Status = "Test")

# inf_prediction <-  observed_contrast_inf |>
#   left_join(predicted_contrast_train, by = c("Comparison", "Feature")) |>
#   group_by(Comparison) |>
#   summarise(
#     # Calculate the MAE
#     MAE = mean(abs(Observed_contrast - Predicted_contrast)),
#     # Calculate the RMSE
#     RMSE = sqrt(mean((Observed_contrast - Predicted_contrast)^2)),
#     # Calculate the R-squared
#     SS_residual = sum((Observed_contrast - Predicted_contrast)^2),  
#     SS_total = sum((Observed_contrast - mean(Observed_contrast))^2), 
#     R2 = 1 - (SS_residual / SS_total),
#     .groups = "drop"
#   ) |>
#   mutate(Status = "Inference")

# # Combine
# metric_prediction <- bind_rows(test_prediction, inf_prediction) |>
#   select(-Comparison) |>
#   distinct() |>
#   left_join(metric_df, by = "Status")

# # Save
# fwrite(metric_prediction, file.path(setup$output_directory, "Metric_contrast.csv"))


# Predicting regression (Baseline differences) ----------------------------------------------------------------------
print("Predicting baseline differences.")
train_coefficients <- train_limma_fit$coefficients
# Prediction
test_predictions <- as_tibble(t(train_coefficients %*% t(test_design)), rownames = "Idx")
inf_predictions <- as_tibble(t(train_coefficients %*% t(inf_design)), rownames = "Idx")
# Observed
test_observations <- as_tibble(t(test_norm$E), rownames = "Idx")
inf_observations <- as_tibble(t(inf_norm$E), rownames = "Idx")
# Meta
test_meta <- as_tibble(test_filtered$obs, rownames = "Idx") |> select(Idx, setup$classification$output_column)
inf_meta <- as_tibble(inf_filtered$obs, rownames = "Idx") |> select(Idx, setup$classification$output_column)
# Merge
test_predictions <- left_join(test_meta, test_predictions, by = "Idx")
inf_predictions <- left_join(inf_meta, inf_predictions, by = "Idx")

test_observations <- left_join(test_meta, test_observations, by = "Idx")
inf_observations <- left_join(inf_meta, inf_observations, by = "Idx")

# Reshape
test_observations <- test_observations |>
  pivot_longer(
    cols = -c(Idx, setup$classification$output_column), 
    names_to = "Feature", 
    values_to = "Observed"
    )

inf_observations <- inf_observations |>
  pivot_longer(
    cols = -c(Idx, setup$classification$output_column), 
    names_to = "Feature", 
    values_to = "Observed"
    )

test_predictions <- test_predictions |>
  pivot_longer(
    cols = -c(Idx, setup$classification$output_column), 
    names_to = "Feature", 
    values_to = "Predicted"
    )

inf_predictions <- inf_predictions |>
  pivot_longer(
    cols = -c(Idx, setup$classification$output_column), 
    names_to = "Feature", 
    values_to = "Predicted"
    )

# Merge (Observed, predicted)
test_obs_pred <- left_join(test_observations, test_predictions, by = c("Idx", setup$classification$output_column, "Feature"))
inf_obs_pred <- left_join(inf_observations, inf_predictions, by = c("Idx", setup$classification$output_column, "Feature"))

# ---- Per gene metric ----
test_gene_error <- test_obs_pred |>
  group_by(!!sym(setup$classification$output_column), Feature) |>
  summarise(
    RMSE = sqrt(mean((Observed - Predicted)^2)),  
    MAE  = mean(abs(Observed - Predicted)),   
    R2 = 1 - sum((Observed - Predicted)^2) / sum((Observed - mean(Observed))^2),   
    .groups = "drop"
  ) |>
  mutate(
    Status = "Test",
    Prediction = "Subset"
  )

inf_gene_error <- inf_obs_pred |>
  group_by(!!sym(setup$classification$output_column), Feature) |>
  summarise(
    RMSE = sqrt(mean((Observed - Predicted)^2)),  
    MAE  = mean(abs(Observed - Predicted)),
    R2 = 1 - sum((Observed - Predicted)^2) / sum((Observed - mean(Observed))^2),  
    .groups = "drop"
  ) |>
  mutate(
    Status = "Inference",
    Prediction = "Ancestry"
  )

# Combine
baseline_per_gene_metric <- bind_rows(test_gene_error, inf_gene_error) |>
  mutate(
    Seed = setup$seed,
    Ancestry = toupper(setup$classification$infer_ancestry),
    n_train_ancestry = train_n,
    n_inf_ancestry = inf_n,
  )

baseline_per_gene_metric <- baseline_per_gene_metric |>
  left_join(condition_count, by = setup$classification$output_column) |>
  rename(n_condition = n) |>
  rename(Condition = !!sym(setup$classification$output_column))

# Save
fwrite(baseline_per_gene_metric, file.path(setup$output_directory, "Baseline_metric_per_gene_dge.csv"))

# ---- Summarized metric ----
test_metric <- test_obs_pred |>
  group_by(!!sym(setup$classification$output_column)) |>
  summarise(
    RMSE = sqrt(mean((Observed - Predicted)^2)),  
    MAE  = mean(abs(Observed - Predicted)),      
    R2   = 1 - (sum((Observed - Predicted)^2) / sum((Observed - mean(Observed))^2)),
    .groups = "drop"
  ) |>
  mutate(
    Status = "Test",
    Prediction = "Subset"
    )

inf_metric <- inf_obs_pred |>
  group_by(!!sym(setup$classification$output_column)) |>
  summarise(
    RMSE = sqrt(mean((Observed - Predicted)^2)),  
    MAE  = mean(abs(Observed - Predicted)),      
    R2   = 1 - (sum((Observed - Predicted)^2) / sum((Observed - mean(Observed))^2)),
    .groups = "drop"
  ) |>
  mutate(
    Status = "Inference",
    Prediction = "Ancestry"
  )

# Combine test_metric and inf_metric
baseline_metric <- bind_rows(test_metric, inf_metric) |>
  mutate(
    Seed = setup$seed,
    Ancestry = toupper(setup$classification$infer_ancestry),
    n_train_ancestry = train_n,
    n_inf_ancestry = inf_n,
  )

baseline_metric <- baseline_metric |>
  left_join(condition_count, by = setup$classification$output_column) |>
  rename(n_condition = n) |>
  rename(Condition = !!sym(setup$classification$output_column))

# Save 
fwrite(baseline_metric, file.path(setup$output_directory, "Baseline_metric_dge.csv"))


print("Switching back to Python.")



# meta = inf_data$obs
# signal = inf_norm$E
# mean_stats = inf_mean_res
# contrast_stats = inf_contrast_res
# goi = goi
# data_output_column = setup$classification$output_column



# mean_stats <- 
#     mean_stats |> 
#     as_tibble() |> 
#     group_by(coef) |> 
#     filter(Gene %in% goi) 

# contrast_stats <- 
#     contrast_stats |> 
#     as_tibble() |> 
#     group_by(coef) |> 
#     filter(Gene %in% goi) 
  
# signal <- 
#     signal |> 
#     as_tibble(rownames = 'Gene') |> 
#     filter(Gene %in% goi) |> 
#     as.data.frame() |> 
#     column_to_rownames('Gene') |> 
#     t() |> 
#     as_tibble(rownames = 'sampleId') 
  
# mean_signal <- 
#     meta |> 
#     as_tibble(rownames = 'sampleId') |> 
#     select(all_of(data_output_column), sampleId) |> 
#     merge(signal, by = 'sampleId') |> 
#     group_by_at(data_output_column) |> 
#     summarise(across(where(is.numeric), mean)) |> 
#     pivot_longer(cols = !all_of(data_output_column),
#                  values_to = 'E',
#                  names_to = 'Gene')
  
#   # Plot stats. results
#   contrast_stats_plot <- ggplot(
#     contrast_stats,
#     aes(
#       x = coef,
#       y = Gene,
#       color = logFC,
#       size = pmin(5, -log10(adj.P.Val))
#     )
#   ) +
#        geom_point() +
#        scale_color_gradient2(
#         low="blue",
#         high="red"
#   ) +
#   ggtitle('Results of comparison') +
#   ylab('Gene') +
#   xlab('Data output column') +
#   theme(axis.text.x = element_text(angle = 90))

#     mean_stats_plot <- ggplot(
#     mean_stats,
#     aes(
#       x = coef,
#       y = Gene,
#       color = logFC,
#       size = pmin(5, -log10(adj.P.Val))
#     )
#   ) +
#        geom_point() +
#        scale_color_gradient2(
#         low="blue",
#         high="red"
#   ) +
#   ggtitle('Results of means model') +
#   ylab('Gene') +
#   xlab('Data output column') +
#   theme(axis.text.x = element_text(angle = 90))

#   # Plot mean signal
#   mean_signal_plot <- ggplot(
#     mean_signal,
#     aes(
#       x = get(data_output_column),
#       y = Gene,
#       fill = E
#     )
#   ) +
#   geom_tile() +
#   scale_fill_gradient2(
#         low="blue",
#         high="red"
#   ) +
#   ggtitle('Mean signal in the data') +
#   ylab('Gene') +
#   xlab('Data output column') +
#   theme(axis.text.x = element_text(angle = 90))

#   # Patchwork
#   validation_plot <- contrast_stats_plot + mean_stats_plot + mean_signal_plot
#   validation_plot

# Metrics to score the model
# Design matrix 
# samples x condition

# Coefficient matrix
# genes x condition

# Metric for goodnes of fit
# y_hat <- fitted(train_limmaFit) # values
# residuals <- train_norm$E - y_hat

# # Sum of square total (SST) based on Europeans
# mean_y <- rowMeans(train_norm$E)  # mean expression for each gene across all samples
# sst <- rowSums((train_norm$E - mean_y)^2) # the total error per gene
# # Residuals sum of squares
# ssr <- rowSums(residuals^2)
# # R2
# r2 <- 1 - (ssr / sst)
# mean_r2 <- mean(r2)

# # RMSE
# rmse <- sqrt(rowMeans(residuals^2))
# mean_rmse <- mean(rmse)

# # Apply European model to Asian
# design_europe <- model.matrix(~ cancer_type_detailed, data = train_data$obs)
# fit_europe <- lmFit(train_norm$E, design_europe)
# fit_europe <- eBayes(fit_europe)

# y_hat_europe <- fitted(fit_europe)
# residuals_europe <- train_norm$E - y_hat_europe

# mean_y_europe <- rowMeans(train_norm$E)
# sst_europe <- rowSums((train_norm$E - mean_y_europe)^2)
# sse_europe <- rowSums(residuals_europe^2)
# r2_europe <- 1 - (sse_europe / sst_europe)
# rmse_europe <- sqrt(rowMeans(residuals_europe^2))
# mean(rmse_europe)

# # Asian
# design_asian <- model.matrix(~ cancer_type_detailed, data = inf_data$obs)
# # Matrix multiplication for the right dimensions
# y_hat_asian <- t(design_asian %*% t(fit_europe$coefficients))
# residuals_asian <- inf_norm$E - y_hat_asian

# mean_y_asian <- rowMeans(inf_norm$E)
# sst_asian <- rowSums((inf_norm$E - mean_y_asian)^2)
# sse_asian <- rowSums(residuals_asian^2)
# r2_asian <- 1 - (sse_asian / sst_asian)
# rmse_asian <- sqrt(rowMeans(residuals_asian^2))
# mean(rmse_asian)


# mean(r2_asian)

