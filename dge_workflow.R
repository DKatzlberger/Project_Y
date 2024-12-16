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
print("Envoking R.")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  yaml_file <- args[1]
  # Check if it's a valid YAML file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
} else {
  # Dev settings if no command-line argument provided
  yaml_file <- "dev_settings.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}

print("Start differential gene expression analysis.")

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
train_data <- adata[adata$obs_names %in% train_idx, ]
test_data <- adata[adata$obs_names %in% test_idx, ]
inf_data <- adata[adata$obs_names %in% inf_idx, ]

# Assertion: Check arrays
# Dimensions
stopifnot(length(train_idx) == length(train_data$obs_names))
# Order
stopifnot(all(colnames(t(train_data$X)) ==
                row.names(train_data$obs[colnames(t(train_data$X)), ])))

# Create design matrix
# 1. Extract the covariate, if it exists and has a meaningful value
# 2. Check if covariate is discrete or continues
covariate <- if ("covariate" %in% names(setup$classification) && 
                 !is.null(setup$classification$covariate) &&
                 setup$classification$covariate != "") {
  setup$classification$covariate
} else {
  NULL
}

# TODO - Check if in each group is at least 2 female/male 
check_gender_balance(train_data)

train_design <- create_design(setup$classification$output_column,
                              meta = train_data$obs,
                              covariate = covariate)                            
test_design <- create_design(setup$classification$output_column,
                             meta = test_data$obs,
                             covariate = covariate)
inf_design <- create_design(setup$classification$output_column,
                            meta = inf_data$obs,
                            covariate = covariate)


# Create contrast matrix (Only needs to be created once)
# Used for hypothesis testing between groups
contrast_matrix <- create_contrast(colnames(train_design))


# Filter genes by expression (Use the train set)
# TODO - will maybe be outsourced to keep same genes across ancestries
keeper_genes <- filterByExpr(t(train_data$X), train_design)
# Subset
train_filtered <- train_data[, keeper_genes]
test_filtered <- test_data[, keeper_genes]
inf_filtered <- inf_data[, keeper_genes]

# Save feature names (For use in ML)
write_yaml(train_filtered$var_names,
           file.path(setup$output_directory, "Features.yml"))

# Assertion: Check array
# Same genes in all sets
stopifnot(all(train_data$var_names == test_data$var_names))
# Dimensions
stopifnot(table(keeper_genes)[2] == dim(train_filtered$X)[2])

# TODO - DGEList object
# TODO - CalcNormFactors

# Normalization (logCPM)
train_norm <- voom(t(train_filtered$X), train_design, plot = FALSE)
test_norm <- voom(t(test_filtered$X), test_design, plot = FALSE)
inf_norm <- voom(t(inf_filtered$X), inf_design, plot = FALSE)

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

# Assertion: Check if all genes have been tested
stopifnot(all(dim(train_mean_res)[1] / length(unique(train_mean_res$coef)) ==
                table(keeper_genes)[2]))

# Save results
fwrite(train_mean_res, file.path(setup$output_directory, "Means_train.csv"))
fwrite(test_mean_res, file.path(setup$output_directory, "Means_test.csv"))
fwrite(inf_mean_res, file.path(setup$output_directory, "Means_inf.csv"))

# Fit contrast
train_limma_fit_contrast <- contrasts.fit(train_limma_fit, contrast_matrix)
test_limma_fit_contrast <- contrasts.fit(test_limma_fit, contrast_matrix)
inf_limma_fit_contrast <- contrasts.fit(inf_limma_fit, contrast_matrix)

# Ebayes
train_limma_fit_contrast <- eBayes(train_limma_fit_contrast)
test_limma_fit_contrast <- eBayes(test_limma_fit_contrast)
inf_limma_fit_contrast <- eBayes(inf_limma_fit_contrast)

# Results contrast
train_contrast_res <- extract_results(train_limma_fit_contrast)
test_contrast_res <- extract_results(test_limma_fit_contrast)
inf_contrast_res <- extract_results(inf_limma_fit_contrast)

# Save results
fwrite(train_contrast_res,
       file.path(setup$output_directory, "Contrast_train.csv"))
fwrite(test_contrast_res,
       file.path(setup$output_directory, "Contrast_test.csv"))
fwrite(inf_contrast_res,
       file.path(setup$output_directory, "Contrast_inf.csv"))

# Visualization
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

# Calculate metric
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

# Correlation
pearson <- head(compute_correlation(data = merged_log_FC, method = "pearson"), 2)
spearman <- head(compute_correlation(data = merged_log_FC, method = "spearman"), 2)

# Make metric dataframe
metric_df <- inner_join(pearson, spearman, by = c("V1", "V2")) |>
  mutate(
    Seed = setup$seed,
    Status = ifelse(V1 == "logFC_test", "Test", "Inference"),
    Ancestry = toupper(setup$classification$infer_ancestry),
    Prediction = ifelse(V1 == "logFC_test", "Subset", "Ancestry"),
    n_ancestry = inf_n
  )

# Save metric Dataframe
fwrite(metric_df, file.path(setup$output_directory, "Metric_dge.csv"))


# Print statement to switch back to python script
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

