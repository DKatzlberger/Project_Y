# Remove start up messages
suppressPackageStartupMessages(
  {
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # DGE workflow and functional analysis
    library(edgeR)
    library(fgsea)
    # Parallelization
    library(parallel)
    library(furrr)
    # Visualization
    library(patchwork)
    library(ggrepel)
    library(ComplexHeatmap)
    library(GetoptLong)
    library(circlize)
    library(ggVennDiagram)
    # Standard libraries
    library(uuid)
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
    library(glue)
  }
)
# Custom functions
source("r_utils.R")
source("figure_themes.R")


# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  YAML_FILE <- args[1]
  # Check if it's a valid YAML file
  is_yaml_file(YAML_FILE)
} else {
  # Dev settings if no command-line argument provided
  cat("Running interactive mode for development. \n")
  # "example_settings_interactions.yaml"
  YAML_FILE <- "example_settings_cross_ancestry_dge.yaml"
  is_yaml_file(YAML_FILE)
}


# Input
# --- Settings
# Default
DEFAULT_FILE  <- "default_settings_cross_ancestry_dge.yaml"
default_setup <- load_default_settings(DEFAULT_FILE)
# User
user_setup   <- yaml.load_file(YAML_FILE)
# Default and user settings (user will overwrite default)
merged_setup <- deep_merge(default_setup, user_setup)
setup        <- merged_setup$result
log          <- merged_setup$log
print_merge_log(log)

# Check required settings
# Required settings
required_settings <- c(
  # Reproducibility
  "seed",
  # Classification
  "output_column", 
  "class_0", 
  "class_1", 
  # Ancestry
  "ancestry_column", 
  "train_ancestry", 
  "infer_ancestry", 
  # Data
  "data_path", 
  "tech", 
  # Features
  "normalization_method",
  # Output
  "output_directory"
)

check_settings(setup, required_settings)
# Add info to settings
setup$date <- format(as.POSIXlt(Sys.time(), tz = "GMT"), "%Y-%m-%d %H:%M:%S") 
setup$id   <- toupper(substr(UUIDgenerate(), 1, 10))

# Create output directory
path_to_save_location <- setup$output_directory
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}
# Save the settings
save_name <- file.path(path_to_save_location, "Settings.yaml")
write_yaml(setup, save_name)

# --- Data
output_column     <- setup$output_column
class_0           <- setup$class_0
class_1           <- setup$class_1
ancestry_column   <- setup$ancestry_column
train_ancestry    <- setup$train_ancestry
infer_ancestry    <- setup$infer_ancestry
covariate         <- setup$covariate

# Load data
data_path  <- setup$data_path
is_h5ad_file(data_path)
adata      <- read_h5ad(data_path)
# Check if columns exist in data
required_columns <- c(output_column, ancestry_column, covariate)
check_columns(adata$obs, required_columns)

# Define classification 
comparison <- c(class_0, class_1)
check_values(adata$obs, output_column, comparison)
adata      <- adata[adata$obs[[output_column]] %in% comparison]
# Factorize
adata$obs[[output_column]] <- factor(adata$obs[[output_column]], levels = comparison)

# Define ancestries
ancestries <- c(train_ancestry, infer_ancestry)
check_values(adata$obs, ancestry_column, ancestries)
adata      <- adata[adata$obs[[ancestry_column]] %in% ancestries]
# Factorize
adata$obs[[ancestry_column]] <- factor(adata$obs[[ancestry_column]], levels = ancestries)

# Cross ancestry analysis
# Seed
seed <- setup$seed
set.seed(seed)
# Message
cat(sprintf("New analysis with id: %s; created: %s \n", setup$id, setup$date))
cat(sprintf("Reproducibility seed: %s \n", seed))
cat(sprintf("Save location: %s \n", path_to_save_location))
cat("-------------------------------------------------------------------- \n")


cat("Number of observations. \n")
# Add to settings
counts <- table(adata$obs[[output_column]], adata$obs[[ancestry_column]])
setup$n_class_0_train_ancestry <- counts[class_0, train_ancestry]
setup$n_class_1_train_ancestry <- counts[class_1, train_ancestry]
setup$n_class_0_infer_ancestry <- counts[class_0, infer_ancestry]
setup$n_class_1_infer_ancestry <- counts[class_1, infer_ancestry]


# Visualize: Sample sizes
save_name     <- file.path(path_to_save_location, "QC_sample_sizes.pdf")
# Plot
p_count       <- plot_output_column_count(adata$obs, ancestry_column, output_column)
p_proportions <- plot_output_column_proportion(adata$obs, ancestry_column, output_column)
# Combine
p <- p_count + p_proportions + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
save_ggplot(p, save_name, width = 6, height = 4)
# Print statement
cat("Check plot: 'QC_sample_sizes.pdf' \n")
cat("-------------------------------------------------------------------- \n")

# --- Model formula
# Define groups to compare
adata$obs["group"] <- factor(
  paste(
    adata$obs[[ancestry_column]], 
    adata$obs[[output_column]], 
    sep = "."
  ), 
  levels = c(
    paste(train_ancestry, class_0, sep = "."),  
    paste(train_ancestry, class_1, sep = "."),  
    paste(infer_ancestry, class_0, sep = "."),    
    paste(infer_ancestry, class_1, sep = ".")     
  )
)

# Marix without covariate
formula <- as.formula("~0 + group")

# --- Filter features 
# Settings
tech            <- setup$tech
data_type       <- setup$data_type
filter_features <- setup$filter_features
# Strength of filter
percentile <- setup$percentile

# Name of QC plot
save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
# Titles
title_before <- ggtitle("Before feature filtering")
title_after  <- ggtitle("After feature filtering")
if (filter_features & tech == "transcriptomics"){

  # Print statement
  cat(sprintf("Filtering features (technology: %s). \n", tech))
  cat(sprintf("By count: Threshold '%s' percentile. \n", percentile))

  # Transform to logCPM
  norm_factors <- calculate_tmm_norm_factors(adata$X)
  cpm_data     <- cpm(adata$X, norm_factors = norm_factors, log = TRUE)

  # Filter by signal/count
  min_counts        <- signal_by_percentile(cpm_data, percentile)
  filtered_features <- filter_by_signal(cpm_data, min_counts)
  # Subset
  filtered_data <- adata[, filtered_features]

  # Print statement
  cat(sprintf("By variance: Threshold '%s' percentile. \n", percentile))

  # Transform to logCPM
  norm_factors <- calculate_tmm_norm_factors(filtered_data$X)
  cpm_data     <- cpm(filtered_data$X, norm_factors = norm_factors, log = TRUE)

  # Filter by variance
  min_variance      <- variance_by_percentile(cpm_data, percentile)
  filtered_features <- filter_by_variance(cpm_data, var_threshold = min_variance)
  # Subset
  filtered_data <- filtered_data[, filtered_features]

  # Visualize: Filtering
  data_before <- log2(adata$X + 0.5)
  data_after  <- log2(filtered_data$X + 0.5)

  # Axis
  x_axis <- paste0("log2(", data_type, " + 0.5)")
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis) + title_before
  p_after  <- plot_mean_variance_trend(data_after, x_axis) + title_after
  # Combine
  p <- p_before / p_after
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (filter_features & tech == "methylation") {

  # Transform to mvalues
  mvals <- beta_to_mvalue(adata$X)
  
  # Filter by variance
  min_variance       <- variance_by_percentile(mvals, percentile)
  filtered_features  <- filter_by_variance(mvals, var_threshold = min_variance)
  # Subset
  filtered_data = adata[, filtered_features]
  
  # Visualize: Filtering
  data_before <- adata$X 
  data_after  <- filtered_data$X

  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis) + title_before
  p_after  <- plot_mean_variance_trend(data_after, x_axis) + title_after
  # Combine
  p <- p_before / p_after
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (filter_features & tech == "proteomics"){

  # No filtering
  filtered_data = adata

  # Visualize: Filtering
  data_before <- adata$X 

  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis) + title_before
  # Save
  save_ggplot(p_before, save_name, width = 6, height = 4)

} else{

  # Print statement
  cat("No filtering of features. \n")

  # No filtering
  filtered_data <- adata

  # Visualize: Filtering
  data_before <- adata$X 
  
  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis)
  # Save
  save_ggplot(p_before, save_name, width = 6, height = 4)

}
# Save feature 
save_name <- file.path(path_to_save_location, "Features.yaml")
write_yaml(filtered_data$var_names, save_name)
# Number of features
setup$n_features  <- ncol(filtered_data)
# Message
if (filter_features){
  cat(sprintf("Number of features after filtering: %s (%s). \n", ncol(filtered_data), ncol(adata)))
  cat("Check plot: 'QC_mean_variance_trend.pdf' \n")
  cat("-------------------------------------------------------------------- \n")
}

# --- Ancestry stratification
cat("Creating stratified subset. \n")
train_adata_ <- filtered_data[filtered_data$obs[[ancestry_column]] %in% train_ancestry]
infer_adata  <- filtered_data[filtered_data$obs[[ancestry_column]] %in% infer_ancestry]

# Stratified subset
strata  <- table(infer_adata$obs[output_column])
indices <- stratified_subset(train_adata_$obs, strata, output_column, seed)
# Check data leakage
stopifnot(!length(intersect(indices$train_idx, indices$test_idx)) > 0)
# Subsetting
train_adata <- train_adata_[indices$train_idx, ]
test_adata  <- train_adata_[indices$test_idx, ]

# Visualize: Stratification
train_meta <- train_adata$obs
test_meta  <- test_adata$obs
infer_meta <- infer_adata$obs
# Add information
train_meta$Set <- "Train"
test_meta$Set  <- "Test"
infer_meta$Set <- "Infer"
# Combine
meta_ <- bind_rows(train_meta, test_meta, infer_meta)
# Plot
p_count       <- plot_output_column_count(meta_, "Set", output_column)
p_proportions <- plot_output_column_proportion(meta_, "Set", output_column)
# Combine
p <- p_count + p_proportions + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
save_name <- file.path(path_to_save_location, "QC_ancestry_stratification.pdf")
save_ggplot(p, save_name, width = 6, height = 4)
# Print statement
cat("Check plot: 'QC_ancestry_stratification.pdf' \n")
cat("-------------------------------------------------------------------- \n")


# --- Workflow: Non-bootstrap (Ground truth)
# Design matrix
cat("Calculating observed logFC. \n")
# Matrix
train_design  <- make_clean_design(train_adata$obs, formula)
test_design   <- make_clean_design(test_adata$obs, formula)
infer_design  <- make_clean_design(infer_adata$obs, formula)

# Normalization/Transformation
# Select normalization method
tech                  <- setup$tech
normalization         <- setup$normalization_method
normalization_method  <- normalization_methods[[tech]][[normalization]]$"function"
values_output_name    <- normalization_methods[[tech]][[normalization]]$"output_name"
# Transpose (rows = Genes, cols = Samples)
train_data_t <- t(train_adata$X)
test_data_t  <- t(test_adata$X)
infer_data_t <- t(infer_adata$X)

# Normalization
train_norm <- normalization_method(train_data_t, train_design)
test_norm  <- normalization_method(test_data_t, test_design)
infer_norm <- normalization_method(infer_data_t, infer_design)
# Extract normalized matrix (used for plotting)
train_norm_matrix <- if (is.list(train_norm) && !is.null(train_norm$E)) train_norm$E else train_norm
test_norm_matrix  <- if (is.list(test_norm) && !is.null(test_norm$E)) test_norm$E else test_norm
infer_norm_matrix <- if (is.list(infer_norm) && !is.null(infer_norm$E)) infer_norm$E else infer_norm

# Visualization: Normalized features
vis_features <- setup$vis_features
if (is.null(vis_features)){
  vis_features <- train_adata$var_names[1:9]
} else{
  vis_features <- vis_features
}

# Means model
train_limma_fit <- lmFit(train_norm, design = train_design)
test_limma_fit  <- lmFit(test_norm, design = test_design)
infer_limma_fit <- lmFit(infer_norm, design = infer_design)
# Ebayes
train_limma_fit <- eBayes(train_limma_fit)
test_limma_fit  <- eBayes(test_limma_fit)
infer_limma_fit <- eBayes(infer_limma_fit)

# Contrast
train_calculations <- paste(colnames(train_design)[1], "-", colnames(train_design)[2])
test_calculations  <- paste(colnames(test_design)[1], "-", colnames(test_design)[2])
infer_calculations <- paste(colnames(infer_design)[1], "-", colnames(infer_design)[2])
# Print statement
cat(sprintf("Formula:       %s     \n", deparse(formula)))
cat(sprintf("Groups:        %s%19s \n", paste(colnames(train_design), collapse = paste0(" ")), "Train/Test"))
cat(sprintf("               %s%10s \n", paste(colnames(infer_design), collapse = paste0(" ")), "Infer"))
cat(sprintf("Calculations:  %s%17s \n", train_calculations, "Train/Test"))
cat(sprintf("               %s%9s  \n", infer_calculations, "Infer "))

# Create contrast matrix
train_contrast_matrix <- makeContrasts(contrasts = train_calculations, levels = train_design)
test_contrast_matrix  <- makeContrasts(contrasts = test_calculations, levels = test_design)
infer_contrast_matrix <- makeContrasts(contrasts = infer_calculations, levels = infer_design)

# Fit contrast
train_limma_fit_contrast <- contrasts.fit(train_limma_fit, train_contrast_matrix)
test_limma_fit_contrast  <- contrasts.fit(test_limma_fit, test_contrast_matrix)
infer_limma_fit_contrast <- contrasts.fit(infer_limma_fit, infer_contrast_matrix)

# Ebayes
train_limma_fit_contrast <- eBayes(train_limma_fit_contrast)
test_limma_fit_contrast  <- eBayes(test_limma_fit_contrast)
infer_limma_fit_contrast <- eBayes(infer_limma_fit_contrast)

# Results contrast
train_contrast_res <- extract_results(train_limma_fit_contrast)
test_contrast_res  <- extract_results(test_limma_fit_contrast)
infer_contrast_res <- extract_results(infer_limma_fit_contrast)

# Some information
# Add status (information on subset)
train_contrast_res$Set <- "Train"
test_contrast_res$Set  <- "Test"
infer_contrast_res$Set <- "Inference"
# Add subset information
train_contrast_res$Subset <- seed
test_contrast_res$Subset  <- seed
infer_contrast_res$Subset <- seed
# Observed
train_contrast_res$Statistic <- "Observed"
test_contrast_res$Statistic  <- "Observed"
infer_contrast_res$Statistic <- "Observed"

# Save
train_save_name <- file.path(path_to_save_location, "LogFC_observed_train.csv")
test_save_name  <- file.path(path_to_save_location, "LogFC_observed_test.csv")
infer_save_name <- file.path(path_to_save_location, "LogFC_observed_infer.csv")

fwrite(train_contrast_res, train_save_name)
fwrite(test_contrast_res, test_save_name)
fwrite(infer_contrast_res, infer_save_name)
cat("Check results: 'LogFC_observed_train.csv' 'LogFC_observed_test.csv' 'LogFC_observed_infer.csv' \n")
cat("-------------------------------------------------------------------- \n")


# --- Workflow: H0 Bootstrap
cat("Calculating bootstrapped logFC. \n")
# Combine test and infer into one sample
H0_matrix <- rbind(test_adata$X, infer_adata$X)
H0_meta   <- bind_rows(test_adata$obs, infer_adata$obs)
# Recompute group column: class_0 vs class_1 (Ancestry independent)
H0_meta$group <- factor(H0_meta[[output_column]], levels = c(class_0, class_1))

# Bootstrapped samples
n_boots <- setup$n_bootstraps
# Function call
test_bootstrap_res <- perform_H0_bootstrap_scratch(
  matrix        = H0_matrix,
  meta          = H0_meta,
  size          = nrow(test_adata),
  group_column  = "group",
  normalization = normalization_method,
  n_iterations  = 100
)

infer_bootstrap_res <- perform_H0_bootstrap_scratch(
  matrix        = H0_matrix,
  meta          = H0_meta,
  size          = nrow(test_adata),
  group_column  = "group",
  normalization = normalization_method,
  n_iterations  = 100
)
cat("-------------------------------------------------------------------- \n")


cat("Calculating test statistics. \n")
# --- Observed statistic
# Convert to data.table
setDT(train_contrast_res)
setDT(test_contrast_res)
setDT(infer_contrast_res)

# Rename columns for clarity
train_obs_stat <- train_contrast_res[, .(Feature, logFC_train = logFC)]
test_obs_stat  <- test_contrast_res[, .(Feature, logFC_test = logFC)]
infer_obs_stat <- infer_contrast_res[, .(Feature, logFC_infer = logFC)]

# Correlation to train set (across genes)
test_obs_cor <- merge(test_obs_stat, train_obs_stat, by = "Feature")[,
  .(Pearson  = cor(logFC_train, logFC_test, method = "pearson"),
  Spearman   = cor(logFC_train, logFC_test, method = "spearman"),
  Set        = "Test")
]

infer_obs_cor <- merge(infer_obs_stat, train_obs_stat, by = "Feature")[,
  .(Pearson = cor(logFC_train, logFC_infer, method = "pearson"),
  Spearman  = cor(logFC_train, logFC_infer, method = "spearman"),
  Set       = "Inference")
]

# Combine and compute difference
obs_cor_diff <- rbind(test_obs_cor, infer_obs_cor)
obs_cor_diff <- dcast(obs_cor_diff, . ~ Set, value.var = c("Pearson", "Spearman"))[, 
  .(DiffPearson = Pearson_Test - Pearson_Inference,
  DiffSpearman  = Spearman_Test - Spearman_Inference,
  Statistic     = "Observed")
]

# Euclidean distances (gene-wise)
test_obs_euclid <- merge(test_obs_stat, train_obs_stat, by = "Feature")[,
  .(Euclidean = sqrt((logFC_train - logFC_test)^2), Set = "Test"),
  by = Feature
]

infer_obs_euclid <- merge(infer_obs_stat, train_obs_stat, by = "Feature")[,
  .(Euclidean = sqrt((logFC_train - logFC_infer)^2), Set = "Inference"),
  by = Feature
]

# Difference in Euclidean distance
obs_euclid_diff <- rbind(test_obs_euclid, infer_obs_euclid)
obs_euclid_diff <- dcast(obs_euclid_diff, Feature ~ Set, value.var = "Euclidean")[, 
  .(Feature,
  Metric    = "Euclidean",
  Difference = Test - Inference,
  Statistic  = "Observed")
]

# Difference in logFC (gene-wise)
obs_logFC_diff <- merge(test_obs_stat, infer_obs_stat, by = "Feature")[, 
  .(Feature,
  DifflogFC = logFC_test - logFC_infer,
  Statistic = "Observed")
]

# --- H0 bootstrap statistic (Null distribution)
# Convert to data.table
setDT(test_bootstrap_res)
setDT(infer_bootstrap_res)
setDT(train_obs_stat)

# Prepare: Rename columns
test_boot_stat  <- test_bootstrap_res[, .(bootstrap, Feature, logFC_test = logFC)]
infer_boot_stat <- infer_bootstrap_res[, .(bootstrap, Feature, logFC_infer = logFC)]

# Correlation to train set
cat("Correlation to train set. \n")
test_boot_cor <- merge(test_boot_stat, train_obs_stat, by = "Feature")[,
  .(Pearson  = cor(logFC_train, logFC_test, method = "pearson"),
    Spearman = cor(logFC_train, logFC_test, method = "spearman"),
    Set      = "Test"), 
  by = bootstrap
]

infer_boot_cor <- merge(infer_boot_stat, train_obs_stat, by = "Feature")[,
  .(Pearson  = cor(logFC_train, logFC_infer, method = "pearson"),
    Spearman = cor(logFC_train, logFC_infer, method = "spearman"),
    Set      = "Inference"), 
  by = bootstrap
]

# Difference in correlation
boot_cor_diff <- rbind(test_boot_cor, infer_boot_cor)
boot_cor_diff <- dcast(boot_cor_diff, bootstrap ~ Set, value.var = c("Pearson", "Spearman"))
boot_cor_diff[, `:=`(
  DiffPearson  = Pearson_Test - Pearson_Inference,
  DiffSpearman = Spearman_Test - Spearman_Inference,
  Statistic    = "H0_bootstrapped"
)]

# Euclidean distances to trian set
cat("Euclidean distance to train set. \n")
test_boot_euclid <- merge(test_boot_stat, train_obs_stat, by = "Feature")[,
  .(Euclidean = sqrt((logFC_train - logFC_test)^2), Set = "Test"), 
  by = .(bootstrap, Feature)
]

infer_boot_euclid <- merge(infer_boot_stat, train_obs_stat, by = "Feature")[,
  .(Euclidean = sqrt((logFC_train - logFC_infer)^2), Set = "Inference"), 
  by = .(bootstrap, Feature)
]

boot_euclid_diff <- rbind(test_boot_euclid, infer_boot_euclid)
boot_euclid_diff <- dcast(boot_euclid_diff, bootstrap + Feature ~ Set, value.var = "Euclidean")
boot_euclid_diff[, `:=`(
  Metric     = "Euclidean",
  Difference = Test - Inference,
  Statistic  = "H0_bootstrapped"
)]

# LogFC difference
cat("LogFC differences between subset and ancestry. \n")
boot_logFC_diff <- merge(test_boot_stat, infer_boot_stat, by = c("bootstrap", "Feature"))
boot_logFC_diff[, `:=`(
  DifflogFC = logFC_test - logFC_infer,
  Statistic = "H0_bootstrapped"
)]

# # --- Visualize: Null distribution
# # Difference in correlation
# boot_cor_diff_combined <- boot_cor_diff |>
#   select(-bootstrap) |>
#   bind_rows(obs_cor_diff) |>
#   rename(
#     Pearson  = DiffPearson,
#     Spearman = DiffSpearman
#   ) |>
#   pivot_longer(
#     cols      = c("Pearson", "Spearman"),
#     names_to  = "Metric",
#     values_to = "Value"
#   ) 

# # Plot
# p <- plot_cor_diff_histogram(boot_cor_diff_combined, x = "Value")
# # Save
# save_name <- file.path(path_to_save_location, "QC_correlation_difference_null_dist.pdf")
# save_ggplot(p, save_name, width = 6, height = 4)

# # Difference in euclidean
# boot_euclid_diff_combined <- boot_euclid_diff |>
#   select(-bootstrap) |>
#   bind_rows(obs_euclid_diff) |>
#   pivot_longer(
#     cols      = c("DiffEuclidean"),
#     names_to  = "Metric",
#     values_to = "Value"
#   ) |>
#   filter(Feature %in% vis_features)

# # Plot
# p <- plot_euclid_diff_histogram(boot_euclid_diff_combined, x = "Value")
# # Save
# save_name <- file.path(path_to_save_location, "QC_euclidean_difference_null_dist.pdf")
# save_ggplot(p, save_name, width = 6, height = 4)


# # Difference in logFC
# boot_logFC_diff_combined <- boot_logFC_diff |>
#   bind_rows(obs_logFC_diff) |>
#   pivot_longer(
#     cols      = c("DifflogFC"),
#     names_to  = "Metric",
#     values_to = "Value"
#   ) |>
#   filter(Feature %in% vis_features)

# # Plot
# p <- plot_logFC_diff_histogram(boot_logFC_diff_combined, x = "Value")
# # Save
# save_name <- file.path(path_to_save_location, "QC_logFC_difference_null_dist.pdf")
# save_ggplot(p, save_name, width = 6, height = 4)
# # Print statement
# cat("Check plot: 'QC_correlation_difference_null_dist.pdf', 'QC_logFC_difference_null_dist.pdf' \n")
# cat("-------------------------------------------------------------------- \n")

# --- Transforming bootstrap distribution to paramtetric normal distribution
source("r_utils.R")
null_dist <- compute_pvalues(boot_euclid_diff, obs_euclid_diff, is_global = FALSE)


# Compute null distribution parameters
null_params <- boot_euclid_diff[, .(
  mean_null = mean(Difference),
  var_null  = var(Difference),
  sd_null   = sqrt(var(Difference)),
  n_boot    = .N
), by = Feature]

# Merge with bootstrapped differences
null_dist <- merge(boot_euclid_diff, null_params, by = "Feature")

# Merge observed values
null_dist <- merge(null_dist, obs_euclid_diff[, .(Feature, Observed = Difference)], by = "Feature")

# ---- Pvalue
null_dist[, `:=`(
  p_parametric = 2 * pnorm(-abs(Observed[1] - mean_null[1]) / sqrt(var_null[1])),
  p_empirical  = (sum(abs(Difference) >= abs(Observed[1])) + 1) / (.N + 1)
), by = Feature]

# Summary statistic
summary_statistic <- unique(null_dist[, .(Feature, Observed, Metric, p_parametric, p_empirical)])
summary_statistic[, `:=`(
  padj_parametric = p.adjust(p_parametric, method = "BH"),
  padj_empirical  = p.adjust(p_empirical, method = "BH")
)]

# Sig features
summary_statistic[p_parametric < 0.05]$Feature

# 5. Filter to desired features
plot_data <-  merge(null_dist, summary_statistic[, .(Feature, padj_parametric, padj_empirical)], by = "Feature", all.x = TRUE)

# 6. Generate density curve per gene
density_data <- plot_data[, {
  mu <- mean_null[1]
  sigma <- sd_null[1]
  x_vals <- seq(mu - 4*sigma, mu + 4*sigma, length.out = 200)
  y_vals <- dnorm(x_vals, mean = mu, sd = sigma)
  .(x = x_vals, y = y_vals)
}, by = Feature]

# 7. Extract one observed value per gene
obs_lines <- unique(plot_data[, .(Feature, Observed)])

plot_data[, density_label := "Bootstrap distribution"]
density_data[, density_label := "Fitted normal distribution"]
obs_lines[, obs_label := "Observed value"]

# 8. Prepare p-value labels (fixed position for all facets)
pvals_annot <- unique(plot_data[, .(Feature, p_parametric, p_empirical)])

pvals_annot[, p_parametric_label := sprintf("Parametric p = %.3g", p_parametric)]
pvals_annot[, p_empirical_label := sprintf("Empirical p = %.3g", p_empirical)]
pvals_annot[, p_label := paste(p_parametric_label, p_empirical_label, sep = "\n")]

p <- ggplot() +
  geom_histogram(
    data = plot_data,
    aes(
      x    = Difference, 
      y    = after_stat(density)
    ),
    bins  = 30,
  ) +
  geom_line(
    data = density_data,
    aes(
      x     = x,
      y     = y, 
      color = density_label
    ),
    linewidth = 0.5 / 1.5
  ) +
  geom_vline(
    data = obs_lines,
    aes(
      xintercept = Observed, 
      linetype   = obs_label, 
      color      = obs_label
    ),
    linewidth = 0.5 / 1.5,
    linetype  = "dashed"
  ) +
  geom_text(
    data = pvals_annot,
    aes(
      x     = -Inf,
      y     = Inf,
      label = p_label
    ),
    hjust = -0.1,  
    vjust = 1.1,  
    size  = 0.5 * 3,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Feature) +
  scale_color_manual(
    name = "Line description",
    values = c(
      "Fitted normal distribution" = "blue",
      "Observed value" = "orange"
    )
  ) +
  labs(
    x       = "Difference",
    y       = "Density",
    caption = "Two-sided pvalue: Is there a difference in either direction"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend()

# 9. Save the plot (make sure path exists)
save_name <- file.path(path_to_save_location, "QC_test.pdf")
ggsave(save_name, p, width = 6, height = 4)


# --- Variance shrinkage
euclid_null_var <- boot_euclid_diff[, .(
  var_null = var(DiffEuclidean),
  sd_null  = sd(DiffEuclidean),   
  n_boot   = .N
), by = Feature]
# Here: Global per genes could be based on libraries like Hallmark
squeezed_var <- squeezeVar(euclid_null_var$var_null, df = euclid_null_var$n_boot - 1)
euclid_null_var[, `:=`(
  var_shrunk = squeezed$var.post,
  sd_shrunk  = sqrt(squeezed$var.post)
)]
# Scaled and original null distribution
boot_euclid_diff <- merge(boot_euclid_diff, euclid_null_var[, .(Feature, sd_shrunk)], by = "Feature")
boot_euclid_diff[, DiffEuclidean_scaled := DiffEuclidean / sd_shrunk]

# Visualize
boot_long <- rbind(obs_euclid_diff, boot_euclid_diff, fill = TRUE)

boot_long <- melt(
  boot_long,
  measure.vars  = c("DiffEuclidean", "DiffEuclidean_scaled"),
  variable.name = "Metric",
  value.name    = "Value"
)

boot_long[, Metric := fcase(
  Metric == "DiffEuclidean", "Original",
  Metric == "DiffEuclidean_scaled", "Scaled"
)]

source("figure_themes.R")
p <- plot_euclid_diff_histogram_scaled(boot_long[Feature %in% vis_features])
# Save
save_name <- file.path(path_to_save_location, "QC_test.pdf")
save_ggplot(p, save_name, width = 6, height = 4)


# --- Empirical p-value:
cat("Calculating empirical p-value. \n")
# Difference in correlation
summary_stats_cor <- bind_rows(
    x = obs_cor_diff, 
    y = boot_cor_diff
  ) |>
  pivot_longer(
    cols      = c(DiffPearson, DiffSpearman),
    names_to  = "Metric",
    values_to = "Difference"
  ) |>
  mutate(
    Metric = case_when(
      Metric == "DiffPearson"  ~ "Pearson",
      Metric == "DiffSpearman" ~ "Spearman"
    )
  ) |>
  group_by(Metric) |>
  summarise(
    obs_diff   = Difference[Statistic == "Observed"],
    boot_diffs = list(Difference[Statistic == "H0_bootstrapped"]),
    .groups    = "drop"
  ) |>
  rowwise() |>
  mutate(
    p_subset    = mean(boot_diffs >= obs_diff),                  # Subset correlation > Ancestry
    p_ancestry  = mean(boot_diffs <= obs_diff),                  # Ancestry correlation > Subset
    p_two_sided = mean(abs(boot_diffs) >= abs(obs_diff))         # Two-sided
  ) |>
  ungroup() |> 
  rename(Difference = obs_diff) |>
  select(-boot_diffs) |>
  mutate(
    padj_subset    = p.adjust(p_subset, method = "BH"),
    padj_ancestry  = p.adjust(p_ancestry, method = "BH"),
    padj_two_sided = p.adjust(p_two_sided, method = "BH")
  ) |>
  mutate(
    Regulation = case_when(
      Difference > 0 ~ "Subset lower distance than ancestry",
      Difference < 0 ~ "Subset higher distance than ancestry",
      TRUE           ~ "No difference"
    )
  )

# Difference in euclidean distance
summary_stats_euclidean <- bind_rows(
    x = obs_euclid_diff, 
    y = boot_euclid_diff
  ) |>
  pivot_longer(
    cols      = c(DiffEuclidean),
    names_to  = "Metric",
    values_to = "Difference"
  ) |>
  group_by(Metric, Feature) |>
  summarise(
    obs_diff   = Difference[Statistic == "Observed"],
    boot_diffs = list(Difference[Statistic == "H0_bootstrapped"]),
    .groups    = "drop"
  ) |>
  rowwise() |>
  mutate(
    p_subset    = mean(boot_diffs <= obs_diff),          # Subset closer to Train than Inference
    p_ancestry  = mean(boot_diffs >= obs_diff),          # Ancestry closer to Train than Inference
    p_two_sided = mean(abs(boot_diffs) >= abs(obs_diff)) # Any large deviation
  ) |>
  ungroup() |> 
  rename(Difference = obs_diff) |>
  select(-boot_diffs) |>
  mutate(
    padj_subset    = p.adjust(p_subset, method = "BH"),
    padj_ancestry  = p.adjust(p_ancestry, method = "BH"),
    padj_two_sided = p.adjust(p_two_sided, method = "BH")
  ) |>
  mutate(
    Metric     = "Euclidean",
    Regulation = case_when(
      Difference > 0 ~ "Subest higher distance than ancestry",
      Difference < 0 ~ "Subset lower distance than ancestry",
      TRUE           ~ "No difference"
    )
  )

# Difference in logFC
summary_stats_logFC <- bind_rows(
    x = obs_logFC_diff, 
    y = boot_logFC_diff
  ) |>
  pivot_longer(
    cols      = c(DifflogFC),
    names_to  = "Metric",
    values_to = "Difference"
  ) |>
  group_by(Metric, Feature) |>
  summarise(
    obs_diff   = Difference[Statistic == "Observed"],
    boot_diffs = list(Difference[Statistic == "H0_bootstrapped"]),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(
    p_subset    = mean(boot_diffs >= obs_diff),                  # Subset upregulated
    p_ancestry  = mean(boot_diffs <= obs_diff),                  # Ancestry upregulated
    p_two_sided = mean(abs(boot_diffs) >= abs(obs_diff))         # Two-sided
  ) |>
  ungroup() |> 
  rename(Difference = obs_diff) |>
  select(-boot_diffs) |>
  mutate(
    padj_subset    = p.adjust(p_subset, method = "BH"),
    padj_ancestry  = p.adjust(p_ancestry, method = "BH"),
    padj_two_sided = p.adjust(p_two_sided, method = "BH")
  ) |>
  mutate(
    Metric     = "Distance",
    Regulation = case_when(
      Difference > 0 ~ "Upregulated in Subset",
      Difference < 0 ~ "Upregulted in Ancestry",
      TRUE           ~ "No difference"
    )
  )

# Save
cor_save_name    <- file.path(path_to_save_location, "Summary_stats_cor.csv")
euclid_save_name <- file.path(path_to_save_location, "Summary_stats_euclid.csv")
logFC_save_name  <- file.path(path_to_save_location, "Summary_stats_logFC.csv")

fwrite(summary_stats_cor, cor_save_name)
fwrite(summary_stats_euclid, euclid_save_name)
fwrite(summary_stats_logFC, logFC_save_name)
cat("Check results: 'Summary_stats_cor.csv' 'Summary_stats_euclid.csv' 'Summary_stats_logFC.csv'. \n")
cat("-------------------------------------------------------------------- \n")

# # --- Visualize: P-value distribution
# # Euclidean distance
# euclidean_p_dist <- summary_stats_euclidean |>
#   select(Feature, Difference, p_two_sided, p_subset, p_ancestry) |>
#   pivot_longer(
#     cols      = starts_with("p"),
#     names_to  = "Type",
#     values_to = "Value"
#   )

# p <- ggplot(
#     data = euclidean_p_dist,
#     aes(
#       x    = Value,
#       fill = factor(floor(Difference))
#     )
#   ) +
#   geom_histogram(
#     binwidth = 0.01
#   ) +
#   facet_wrap(
#     ~ Type
#   ) +
#   scale_fill_discrete(
#     name = "Observed statisic"
#   ) + 
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme_small_legend()

# # Save
# save_name <- file.path(path_to_save_location, "QC_euclidean_p_value_dist.pdf")
# save_ggplot(p, save_name, width = 6, height = 4)




# --- Visualize
if (setup$visual_val){

  # Print statement
  cat("Visualizing results. \n")
  cat("-------------------------------------------------------------------- \n")
  # Create output directory
  path_to_save_location <- file.path(setup$output_directory, "Visual_val")
  if (!dir.exists(path_to_save_location)) {
    dir.create(path_to_save_location, recursive = TRUE)
  }

  # --- Correlation
  plot_cor <- bind_rows(
      x = test_obs_cor,
      y = infer_obs_cor
    ) |>
    mutate(
      Prediction = case_when(
        Set == "Test"      ~ "Subset",
        Set == "Inference" ~ "Ancestry",
      )
    ) |>
    pivot_longer(
      cols      = c("Pearson", "Spearman"),
      names_to  = "Metric",
      values_to = "Value"
    ) |>
    inner_join(
      summary_stats_cor,
      by = "Metric"
    )
  
  # --- Euclidean
  # Volcano plot
  p <- plot_euclid_volcano(
    data = summary_stats_euclidean,
    x           = "Difference",
    y           = "p_two_sided",
    arrow_left  = "Subset closer",
    arrow_right = "Ancestry closer",
    thr         = setup$logFC_thr,
    x_label     = expression(Delta*"Euclidean distance"),
    caption     = expression(Delta*"Euclidean distance:" ~ sqrt((logFC[Train] - logFC[Test])^2) - sqrt((logFC[Train] - logFC[Inference])^2))
  )
  # Save
  save_name <- file.path(path_to_save_location, "Euclidean_distance_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Significant features
  sig   <- filter(summary_stats_euclidean, p_two_sided < 0.05 & abs(Difference) > setup$logFC_thr)
  sig   <- arrange(sig, p_two_sided) |> slice_head(n = 10)
  sig_features <- sig$Feature

  # Visualize: Sig features
  if (length(sig_features) != 0){

    # Features with normalised expression
    train_zscore <- get_expression(
      matrix   = train_norm_matrix,
      meta     = train_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Train")

    test_zscore <- get_expression(
      matrix   = test_norm_matrix,
      meta     = test_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Test")

    infer_zscore <- get_expression(
      matrix   = infer_norm_matrix,
      meta     = infer_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Inference")

    # Combine
    zscore <- bind_rows(train_zscore, test_zscore, infer_zscore)
    zscore <- mutate(zscore, Set = factor(Set, levels = c("Train", "Test", "Inference")))
    # Plot
    p <- plot_euclid_boxplot(
      data = zscore,
      x       = "Set",
      y       = "zscore",
      fill    = output_column,
      x_label = "Set"
    )
    # Save
    save_name <- file.path(path_to_save_location, "Euclidean_distance_boxplot.pdf")
    save_ggplot(p, save_name, width = 6, height = 4)
  }

  # --- LogFC
  # Volcano plot
  p <- plot_euclid_volcano(
    data        = summary_stats_logFC,
    x           = "Difference",
    y           = "p_two_sided",
    arrow_left  = "Ancestry upregulated",
    arrow_right = "Subset upregulated",
    thr         = setup$logFC_thr,
    x_label     = expression(Delta*"logFC"),
    caption     = expression(Delta*"logFC:"~logFC[Test] - logFC[Inference])
  )
  # Save
  save_name <- file.path(path_to_save_location, "Interaction_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Significant features
  sig          <- filter(summary_stats_logFC, p_two_sided < 0.05 & abs(Difference) > setup$logFC_thr)
  sig          <- arrange(sig, p_two_sided) |> slice_head(n = 10)
  sig_features <- sig$Feature

  # Visualize: Sig features
  if (length(sig_features) != 0){

    # Features with normalised expression
    train_zscore <- get_expression(
      matrix   = train_norm_matrix,
      meta     = train_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Train")

    test_zscore <- get_expression(
      matrix   = test_norm_matrix,
      meta     = test_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Test")

    infer_zscore <- get_expression(
      matrix   = infer_norm_matrix,
      meta     = infer_adata$obs,
      features = sig_features
    ) |>
    mutate(Set = "Inference")

    # Combine
    zscore <- bind_rows(train_zscore, test_zscore, infer_zscore)
    zscore <- mutate(zscore, Set = factor(Set, levels = c("Train", "Test", "Inference")))
    # Plot
    p <- plot_euclid_boxplot(
      data = zscore,
      x       = "Set",
      y       = "zscore",
      fill    = output_column,
      x_label = "Set"
    )
    # Save
    save_name <- file.path(path_to_save_location, "Interaction_boxplot.pdf")
    save_ggplot(p, save_name, width = 6, height = 4)
  }
}









