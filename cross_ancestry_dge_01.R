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
  p_before <- plot_mean_variance_trend(data_before, x_axis)
  p_after  <- plot_mean_variance_trend(data_after, x_axis)
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
  p_before <- plot_mean_variance_trend(data_before, x_axis) 
  p_after  <- plot_mean_variance_trend(data_after, x_axis)
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
  p_before <- plot_mean_variance_trend(data_before, x_axis) 
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
train_adata <- filtered_data[filtered_data$obs[[ancestry_column]] %in% train_ancestry]
infer_adata <- filtered_data[filtered_data$obs[[ancestry_column]] %in% infer_ancestry]
# Stratified subset
strata  <- table(infer_adata$obs[output_column])
indices <- stratified_subset(train_adata$obs, strata, output_column, seed)
# Check data leakage
stopifnot(!length(intersect(indices$train_idx, indices$test_idx)) > 0)
# Subsetting
train_adata <- adata[indices$train_idx, ]
test_adata  <- adata[indices$test_idx, ]

# Visualize: Stratification
train_meta <- train_adata$obs
test_meta  <- test_adata$obs
infer_meta <- infer_adata$obs
# Add information
train_meta$Status <- "Train"
test_meta$Status  <- "Test"
infer_meta$Status <- "Infer"
# Combine
meta_ <- bind_rows(train_meta, test_meta, infer_meta)
# Plot
save_name     <- file.path(path_to_save_location, "QC_ancestry_stratification.pdf")
p_count       <- plot_output_column_count(meta_, "Status", output_column)
p_proportions <- plot_output_column_proportion(meta_, "Status", output_column)
# Combine
p <- p_count + p_proportions + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
save_ggplot(p, save_name, width = 6, height = 4)
# Print statement
cat("Check plot: 'QC_ancestry_stratification.pdf' \n")
cat("-------------------------------------------------------------------- \n")


# Workflow: Train


# --- Design matrix
cat("Train workflow. \n")
# Matrix
train_design <- model.matrix(formula, data = train_adata$obs)
# Human and machine readable terms
colnames(train_design) <- gsub("group", "", colnames(train_design))
colnames(train_design) <- gsub("-", "_", colnames(train_design))
# Print statement
cat(sprintf("Formula:       %s\n", deparse(formula)))
cat(sprintf("Groups:        %s\n", paste(colnames(train_design), collapse = paste0(" "))))

# --- Normalization/Transformation
# Select normalization method
tech                  <- setup$tech
normalization         <- setup$normalization_method
normalization_method  <- normalization_methods[[tech]][[normalization]]$"function"
values_output_name    <- normalization_methods[[tech]][[normalization]]$"output_name"
# Transpose (rows = Genes, cols = Samples)
train_data_t <- t(train_adata$X)
# Normalization
train_norm  <- normalization_method(train_data_t, train_design)
# Extract normalized matrix (used for plotting)
train_norm_matrix <- if (is.list(train_norm) && !is.null(train_norm$E)) train_norm$E else train_norm

# --- Means model
# Fit the model (means model)
limma_fit      <- lmFit(train_norm, design = train_design)
limma_fit      <- eBayes(limma_fit)
train_mean_res <- extract_results(limma_fit)
# Some information
train_mean_res$Status <- "Train"
train_mean_res$Seed   <- seed

# Save
# save_name <- file.path(path_to_save_location, "Limma_means_train.csv")
# fwrite(train_mean_res, save_name)
# cat("Check results: 'Limma_means_train.csv' \n")

# --- Contrast
# Calculations
cols                  <- colnames(train_design)
contrast_calculations <- glue("{cols[1]} - {cols[2]}")
# Print statement
cat(sprintf("Calculations:  %s \n", contrast_calculations))

# Create contrast matrix
contrast_matrix <- makeContrasts(
  contrasts = contrast_calculations,
  levels    = train_design
)

# Fit contrast
limma_fit_contrast <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit_contrast <- eBayes(limma_fit_contrast)
train_contrast_res <- extract_results(limma_fit_contrast)
# Some information
train_contrast_res$Status <- "Train"
train_contrast_res$Subset <- seed

# Save
save_name <- file.path(path_to_save_location, "Limma_contrast_train.csv")
fwrite(train_contrast_res, save_name)
cat("Check results: 'Limma_contrast_train.csv'. \n")
cat("-------------------------------------------------------------------- \n")

# Workflow: Bootstrap
cat("Test workflow. \n")
test_contrast_res <- perform_bootstrap_parallel(
  expression_matrix        = test_adata$X, 
  metadata                 = test_adata$obs, 
  formula                  = formula, 
  normalization_method     = normalization_method, 
  n_iterations             = 10
)
# Addition information
test_contrast_res$Status <- "Train"
test_contrast_res$Subset <- seed
# Save
save_name <- file.path(path_to_save_location, "Limma_contrast_test.csv")
fwrite(test_contrast_res, save_name)
cat("Check results: 'Limma_contrast_test.csv'. \n")
cat("-------------------------------------------------------------------- \n")

cat("Infer workflow. \n")
infer_contrast_res <- perform_bootstrap_parallel(
  expression_matrix    = infer_adata$X, 
  metadata             = infer_adata$obs,
  formula              = formula, 
  normalization_method = normalization_method, 
  n_iterations         = 10
)
# Addition information
infer_contrast_res$Status <- "Infer"
infer_contrast_res$Subset <- seed
# Save
save_name <- file.path(path_to_save_location, "Limma_contrast_inf.csv")
fwrite(infer_contrast_res, save_name)
cat("Check results: 'Limma_contrast_inf.csv'. \n")
cat("-------------------------------------------------------------------- \n")

# Save the settings
save_name <- file.path(path_to_save_location, "Settings.yaml")
write_yaml(setup, save_name)

# --- Integrating summary statistic (logFC)
# Local (Feature wise)
train_stats <- select(train_contrast_res, Feature, logFC)
test_stats  <- select(test_contrast_res, bootstrap, Feature, logFC)
infer_stats <- select(infer_contrast_res, bootstrap, Feature, logFC)

# Rename for clarity
train_stats <- rename(train_stats, logFC_train = logFC)
test_stats  <- rename(test_stats, logFC_boot = logFC)
infer_stats <- rename(infer_stats, logFC_boot = logFC)

# Compute difference: test_diff
test_diff <- test_stats |>
  inner_join(train_stats, by = "Feature") |>
  mutate(
    Error     = logFC_boot - logFC_train,
    AbsError = abs(logFC_boot - logFC_train),
    Status    = "Test"
    )

# Compute difference: infer_diff
infer_diff <- infer_stats |>
  inner_join(train_stats, by = "Feature") |>
  mutate(
    Error    = logFC_boot - logFC_train,
    AbsError = abs(logFC_boot - logFC_train),
    Status   = "Infer"
    )

# Combine the two datasets
summary_statisitc <- bind_rows(test_diff, infer_diff)

error_diff <- summary_statisitc %>%
  select(Feature, bootstrap, Status, AbsError) %>%
  pivot_wider(
    names_from = Status,
    values_from = AbsError,
    names_prefix = "AbsError_"
  ) %>%
  mutate(
    AbsError_Diff = AbsError_Infer - AbsError_Test
  )
