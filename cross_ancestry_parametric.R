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
source("cross_ancestry_utils.R")
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


# --- Select normalization method
tech                  <- setup$tech
normalization         <- setup$normalization_method
normalization_method  <- normalization_methods[[tech]][[normalization]]$"function"
values_output_name    <- normalization_methods[[tech]][[normalization]]$"output_name"

# --- Workflow: Non-bootstrap (Ground truth)
train_obs <- calculate_logFC(
  matrix        = train_adata$X,
  meta          = train_adata$obs,
  group_column  = output_column,
  normalization = normalization_method
)
test_obs <- calculate_logFC(
  matrix        = test_adata$X,
  meta          = test_adata$obs,
  group_column  = output_column,
  normalization = normalization_method
)
infer_obs <- calculate_logFC(
  matrix        = infer_adata$X,
  meta          = infer_adata$obs,
  group_column  = output_column,
  normalization = normalization_method
)

# --- Observed statistic
# Prepare data
train_obs <- train_obs[, .(Feature, logFC_train = logFC)]
test_obs  <- test_obs[, .(Feature, logFC_test = logFC)]
infer_obs <- infer_obs[, .(Feature, logFC_infer = logFC)]

obs_logFC_diff    <- compute_statistic("logFC", train_obs, test_obs, infer_obs)
obs_euclid_diff   <- compute_statistic("Euclidean", train_obs, test_obs, infer_obs)
obs_pearson_list  <- compute_statistic("Pearson", train_obs, test_obs, infer_obs)
obs_spearman_list <- compute_statistic("Spearman", train_obs, test_obs, infer_obs)

# Obserced correlation values for plotting
obs_pearson_diff <- obs_pearson_list$correlation_difference
obs_pearson_raw  <- obs_pearson_list$raw_correlation

obs_spearman_diff <- obs_spearman_list$correlation_difference
obs_spearman_raw  <- obs_spearman_list$raw_correlation


# Visualization: Normalized features
vis_features <- setup$vis_features
if (is.null(vis_features)){
  vis_features <- train_adata$var_names[1:9]
} else{
  vis_features <- vis_features
}

# --- Workflow: H0 Bootstrap
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
  n_iterations  = 10
)

infer_bootstrap_res <- perform_H0_bootstrap_scratch(
  matrix        = H0_matrix,
  meta          = H0_meta,
  size          = nrow(test_adata),
  group_column  = "group",
  normalization = normalization_method,
  n_iterations  = 10
)

# --- H0 bootstrap statistic (Null distribution)
# Prepare data
test_boot  <- test_bootstrap_res[, .(bootstrap, Feature, logFC_test = logFC)]
infer_boot <- infer_bootstrap_res[, .(bootstrap, Feature, logFC_infer = logFC)]
# Expand observed train for bootstrap IDs
boot_ids <- unique(test_boot$bootstrap)
train_boot <- train_obs[
  , .(bootstrap = boot_ids), by = .(Feature, logFC_train)
]

boot_logFC_diff    <- compute_statistic("logFC",    train_boot, test_boot, infer_boot)
boot_euclid_diff   <- compute_statistic("Euclidean", train_boot, test_boot, infer_boot)
boot_pearson_list  <- compute_statistic("Pearson",  train_boot, test_boot, infer_boot)
boot_spearman_list <- compute_statistic("Spearman", train_boot, test_boot, infer_boot)

# Obserced correlation values for plotting
boot_pearson_diff  <- boot_pearson_list$correlation_difference
boot_spearman_diff <- boot_spearman_list$correlation_difference


# --- Pvalue
pvals_euclid   <- compute_pvalue(boot_euclid_diff, obs_euclid_diff, is_global = FALSE)
pvals_logfc    <- compute_pvalue(boot_logFC_diff, obs_logFC_diff, is_global = FALSE)
pvals_pearson  <- compute_pvalue(boot_pearson_diff, obs_pearson_diff, is_global = TRUE)
pvals_spearman <- compute_pvalue(boot_spearman_diff, obs_spearman_diff, is_global = TRUE)

# --- QC null distribution
p <- plot_null_distribution(pvals_euclid[Feature %in% vis_features])
# Save
save_name <- file.path(path_to_save_location, "QC_null_dist_euclid.pdf")
ggsave(save_name, p, width = 6, height = 4)

p <- plot_null_distribution(pvals_logfc[Feature %in% vis_features])
# Save
save_name <- file.path(path_to_save_location, "QC_null_dist_logFC.pdf")
ggsave(save_name, p, width = 6, height = 4)

p <- plot_null_distribution(pvals_pearson)
# Save
save_name <- file.path(path_to_save_location, "QC_null_dist_pearson.pdf")
ggsave(save_name, p, width = 6, height = 4)

p <- plot_null_distribution(pvals_spearman)
# Save
save_name <- file.path(path_to_save_location, "QC_null_dist_spearman.pdf")
ggsave(save_name, p, width = 6, height = 4)

# --- Save summary statistics
important_cols <- c("Feature", "Observed", "Metric", "p_parametric", "p_empirical", "padj_parametric", "padj_empirical")
summary_euclid   <- unique(pvals_euclid[, ..important_cols])
summary_logfc    <- unique(pvals_logfc[, ..important_cols])
summary_pearson  <- unique(pvals_pearson[, ..important_cols])
summary_spearman <- unique(pvals_spearman[, ..important_cols])

save_name_euclid   <- file.path(path_to_save_location, "Summary_euclid.csv")
save_name_logfc    <- file.path(path_to_save_location, "Summary_logFC.csv")
save_name_pearson  <- file.path(path_to_save_location, "Summary_pearson.csv")
save_name_spearman <- file.path(path_to_save_location, "Summary_spearman.csv")

fwrite(summary_euclid, save_name_euclid)
fwrite(summary_euclid, save_name_logfc)
fwrite(summary_euclid, save_name_pearson)
fwrite(summary_euclid, save_name_spearman)

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

  # Interactions
  p <- plot_cross_ancestry_volcano(
      data        = summary_logfc,
      x           = "Observed",
      y           = "padj_parametric",
      arrow_left  = "Inference upregulated",
      arrow_right = "Test upregulated",
      thr         = setup$logFC_thr,
      features    = summary_logfc[Feature %in% vis_features]$Feature,
      x_label     = expression(Delta*"logFC"),
      caption     = expression(Delta*"logFC:"~logFC[Test] - logFC[Inference])
    ) +
    ggtitle("Interactions")
  # Save
  save_name <- file.path(path_to_save_location, "Interaction_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Euclidean
  p <- plot_cross_ancestry_volcano(
    data        = summary_euclid,
    x           = "Observed",
    y           = "padj_parametric",
    arrow_left  = "Test closer to train",
    arrow_right = "Subset closer to train",
    thr         = setup$logFC_thr,
    features    = summary_logfc[Feature %in% vis_features]$Feature,
    x_label     = expression(Delta*"Euclidean"),
    caption     = expression(Delta*"Euclidean distance:" ~ sqrt((logFC[Train] - logFC[Test])^2) - sqrt((logFC[Train] - logFC[Inference])^2))
  )
  # Save
  save_name <- file.path(path_to_save_location, "Euclidean_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Correlation
  summary_cor <- rbind(summary_pearson, summary_spearman)
  raw_cor     <- rbind(obs_pearson_raw, obs_spearman_raw)

  # Merge raw and difference
  raw_cor_wide       <- dcast(raw_cor, Metric ~ Set, value.var = "Correlation")
  summary_cor_merged <- merge(summary_cor, raw_cor_wide, by = "Metric")

  p <- plot_correlation_comparison(summary_cor_merged)
  # Save
  save_name <- file.path(path_to_save_location, "Correlation_bar.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

}



# --- Sanity check (interactions with limma)
# Create output directory
path_to_save_location <- file.path(setup$output_directory, "Sanity_check")
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Design
interaction_matrix <- rbind(test_adata$X, infer_adata$X)
interaction_meta   <- bind_rows(test_adata$obs, infer_adata$obs)
# Recompute group column
interaction_meta$group <- factor(
  paste(
    limma_meta[[ancestry_column]], 
    limma_meta[[output_column]], 
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
interaction_design <- model.matrix(formula, data = interaction_meta)
# Change coefficient names
colnames(interaction_design) <- gsub("group", "", colnames(interaction_design))
colnames(interaction_design) <- gsub("-", "_", colnames(interaction_design))

# --- Normalization
# Settings
tech                 <- setup$tech
normalization_setup  <- setup$normalization_method
normalization_method <- normalization_methods[[tech]][[normalization_setup]]$"function"
values_output_name   <- normalization_methods[[tech]][[normalization_setup]]$"output_name"

# Transpose (rows = Genes, cols = Samples)
data_t    <- t(interaction_matrix)
data_norm <- normalization_method(data_t, interaction_design)
# Extract normalized matrix (used for plotting)
data_norm_matrix <- if (is.list(data_norm) && !is.null(data_norm$E)) data_norm$E else data_norm

# Visualize: Normalization
title_before <- ggtitle("Before normaliaztion")
title_after  <- ggtitle("After normaliaztion")
# Density per sample
p_before <- plot_density_of_samples(interaction_matrix, x_axis_label = setup$data_type) + title_before
p_after  <- plot_density_of_samples(t(data_norm_matrix), x_axis_label = values_output_name) + title_after
# Combine
p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
p_name <- file.path(path_to_save_location, "QC_density_normalized_values.pdf")
save_ggplot(p, p_name, width = 3, height = 3)

# Q-Q plots per gene
p_before <- plot_qq_of_genes(interaction_matrix, n_features = 5) + title_before
p_after  <- plot_qq_of_genes(t(data_norm_matrix), n_features = 5) + title_after
# Combine
p <- p_before / p_after 
# Save
save_name <- file.path(path_to_save_location, "QC_qq_normalized_values.pdf")
save_ggplot(p, save_name, width = 6, height = 4)

# Print statement 
cat("Check plot: 'QC_density_normalized_values.pdf' and 'QC_qq_normalized_values.pdf' \n")
cat("-------------------------------------------------------------------- \n")


# --- Differential gene expression
cat("Differential gene expression analysis. \n")
# --- Means model
# Fit the model (means model)
limma_fit      <- lmFit(data_norm, design = interaction_design)
limma_fit      <- eBayes(limma_fit)
mean_model_res <- extract_results(limma_fit)

# Save
save_name <- file.path(path_to_save_location, "Limma_means.csv")
fwrite(mean_model_res, save_name)
cat("Fit means model. \n")
cat("Check results: 'Limma_means.csv'. \n")

# --- Contrast 
# Terms
contrast_terms <- list(
  baseline_1      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}.Baseline"),
  baseline_2      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[2]}.Baseline"),
  relationship_1  = glue("{train_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relationship"),
  relationship_2  = glue("{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relatioship"),
  interaction     = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Interaction")
)

# Calculations
cols <- colnames(interaction_design)
contrast_calculations <- list(
  baseline_1      = glue("{cols[1]} - {cols[3]}"),
  baseline_2      = glue("{cols[2]} - {cols[4]}"),
  relationship_1  = glue("{cols[1]} - {cols[2]}"),
  relationship_2  = glue("{cols[3]} - {cols[4]}"),
  interaction     = glue("({cols[1]} - {cols[2]}) - ({cols[3]} - {cols[4]})")
)
# Create contrast matrix
contrast_matrix <- makeContrasts(
  contrasts = contrast_calculations,
  levels    = interaction_design
)
colnames(contrast_matrix) <- contrast_terms

# Fit contrast
limma_fit_contrast <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit_contrast <- eBayes(limma_fit_contrast)
contrast_res       <- extract_results(limma_fit_contrast)

# Save
save_name <- file.path(path_to_save_location, "Limma_contrast.csv")
fwrite(contrast_res, save_name)
cat("Check results: 'Limma_contrast.csv'. \n")
cat("-------------------------------------------------------------------- \n")

# Significant features
interaction_term <- contrast_terms$interaction
interaction      <- filter(contrast_res, coef %in% interaction_term)
# Filter
sig_interaction <- filter(interaction, adj.P.Val < 0.05 & abs(logFC) > setup$logFC_thr)
sig_features    <- sig_interaction$Feature
# Save
save_name <- file.path(path_to_save_location, "Sig_features.yaml")
write_yaml(sig_features, save_name)

# --- Visualize interactions
new_c <- c(ancestry_column, output_column, "comp_level")
interaction_term <- contrast_terms$interaction
interaction      <- filter(contrast_res, coef %in% interaction_term)
interaction      <- separate(interaction, coef, into = new_c, sep = "\\.", remove = FALSE)
# Add information
interaction[[ancestry_column]] <- str_replace_all(interaction[[ancestry_column]], "_", " ")
interaction[[output_column]]   <- str_replace_all(interaction[[output_column]], "_", " ")

# Visualize: Interaction
# Volcano plot
p <- plot_cross_ancestry_volcano(
    data        = interaction,
    x           = "logFC",
    y           = "adj.P.Val",
    arrow_left  = "Inference upregulated",
    arrow_right = "Test upregulated",
    thr         = setup$logFC_thr,
    features    = summary_logfc[Feature %in% vis_features]$Feature,
    x_label     = "logFC",
    caption     = paste("Limma interaction logFC:", unique(interaction$coef))
  ) +
  ggtitle("Interactions")
# Save
save_name <- file.path(path_to_save_location, "Interaction_volcano.pdf")
save_ggplot(p, save_name, width = 6, height = 4)

# Visualize: Sig interaction
if (length(sig_features) != 0){
  
  # Number of features to visualize
  if (length(sig_features) > 10){
    n_sig_features <- 10
  } else{
    n_sig_features <- length(sig_features)
  }

  exp_list       <- list()
  for(feature in sig_features[1:n_sig_features]){
    # Meta data for each feature
    exp_list[[feature]] <- adata$obs |>
      mutate(zscore = scale(data_norm_matrix[feature,])) |>
      rownames_to_column("idx") |>
      remove_rownames()
  }
  # Features with normalised expression
  exp_zscore <- bind_rows(exp_list, .id = "Feature")

  # Boxplot
  p <- interactions_boxplot(
    exp_zscore, 
    x          = ancestry_column, 
    fill       = output_column,
    point_size = 0.8
    )
  # Save
  save_name <- file.path(path_to_save_location, "Interaction_boxplot.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Heatmap
  interaction_heatmap(
    exp_zscore,
    output_column,
    ancestry_column,
    path_to_save_location
  )
}




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









