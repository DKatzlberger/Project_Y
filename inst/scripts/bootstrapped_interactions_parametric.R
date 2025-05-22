#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────
# Cross-Ancestry Parametric Main Script (Development Mode)
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(anndata)
  library(yaml)
  library(glue)
  library(uuid)
  library(reticulate)
  library(parallel)
  library(purrr)
  library(furrr)
  library(limma)
  library(patchwork)

  use_condaenv(
    condaenv = "/opt/conda/envs/ancestry/bin/python",
    required = TRUE
  )
})

# ───────────────────────────────────────────────────────────────
# Load package function files (Development Mode)
# ───────────────────────────────────────────────────────────────

function_dir <- "R"
function_files <- list.files(
  path       = function_dir,
  full.names = TRUE,
  pattern    = "\\.R$"
)

suppressMessages(
  invisible(lapply(function_files, source))
)

if (length(function_files) == 0) {
  cat("⚠ No R function files found in 'R/' directory.\n")
} else {
  cat("✔ Loaded function files from R/: \n")
  cat(paste0("• ", basename(function_files), collapse = "\n"), "\n")
}


# ───────────────────────────────────────────────────────────────
# YAML settings file (from CLI or fallback)
# ───────────────────────────────────────────────────────────────

YAML_FILE <- yaml_argument_loader(
  default_yaml = "inst/extdata/example_settings_cross_ancestry_parametric.yaml"
)

# ───────────────────────────────────────────────────────────────
# Default settings file
# ───────────────────────────────────────────────────────────────

DEFAULT_FILE <- get_default_settings_path(
  file_name = "default_settings_cross_ancestry_parametric.yaml",
  base_dir  = file.path("inst", "config")
)

# ───────────────────────────────────────────────────────────────
# Merge and load
# ───────────────────────────────────────────────────────────────

setup <- load_and_merge_settings(
  default_path = DEFAULT_FILE,
  user_path    = YAML_FILE
)

# ───────────────────────────────────────────────────────────────
# Check required settings
# ───────────────────────────────────────────────────────────────

required_settings <- c(
  "seed",
  "output_column", 
  "class_0", 
  "class_1", 
  "ancestry_column", 
  "train_ancestry", 
  "infer_ancestry", 
  "data_path", 
  "tech", 
  "normalization_method",
  "output_directory"
)

check_required_settings(
  settings      = setup,
  required_keys = required_settings
)

# ───────────────────────────────────────────────────────────────
# Create output directory and save finalized settings
# ───────────────────────────────────────────────────────────────

path_to_save_location <- setup$output_directory
if (!dir.exists(path_to_save_location)) {
  dir.create(
    path      = path_to_save_location,
    recursive = TRUE
  )
}

save_settings(
  settings   = setup,
  file_path  = file.path(setup$output_directory, "Settings.yaml")
)

# ───────────────────────────────────────────────────────────────
# Load and validate input data
# ───────────────────────────────────────────────────────────────

data_path       <- setup$data_path
output_column   <- setup$output_column
ancestry_column <- setup$ancestry_column

is_h5ad_file(
  file_path = data_path
)

adata <- read_h5ad(
  filename = data_path
)

required_columns <- c(
  output_column,
  ancestry_column
)

check_columns_exist(
  df            = adata$obs,
  required_cols = required_columns
)

# ───────────────────────────────────────────────────────────────
# Define and validate classification labels
# ───────────────────────────────────────────────────────────────

class_0 <- setup$class_0
class_1 <- setup$class_1

comparison <- c(
  class_0,
  class_1
)

check_values_exist(
  df       = adata$obs,
  column   = output_column,
  expected = comparison
)

adata <- adata[
  adata$obs[[output_column]] %in% comparison,
]

adata$obs[[output_column]] <- factor(
  x      = adata$obs[[output_column]],
  levels = comparison
)

# ───────────────────────────────────────────────────────────────
# Define and validate ancestry groups
# ───────────────────────────────────────────────────────────────

train_ancestry  <- setup$train_ancestry
infer_ancestry  <- setup$infer_ancestry

ancestries <- c(
  train_ancestry,
  infer_ancestry
)

check_values_exist(
  df       = adata$obs,
  column   = ancestry_column,
  expected = ancestries
)

adata <- adata[
  adata$obs[[ancestry_column]] %in% ancestries,
]

adata$obs[[ancestry_column]] <- factor(
  x      = adata$obs[[ancestry_column]],
  levels = ancestries
)

# ───────────────────────────────────────────────────────────────
# Set seed
# ───────────────────────────────────────────────────────────────

seed <- setup$seed
set.seed(seed)

# ───────────────────────────────────────────────────────────────
# Visualize: Sample sizes
# ───────────────────────────────────────────────────────────────

# Plot
p_count <- plot_output_column_count(
  data        = adata$obs,
  x           = ancestry_column,
  fill        = output_column,
  point_size  = 0.5
)

p_proportions <- plot_output_column_proportion(
  data        = adata$obs,
  x           = ancestry_column,
  fill        = output_column,
  point_size  = 0.5
)

p <- p_count + p_proportions + patchwork::plot_layout(guides = "collect") & 
  ggplot2::theme(legend.position = "bottom")

save_name <- file.path(path_to_save_location, "QC_sample_sizes.pdf")

size <- estimate_plot_size(p)

save_ggplot(
  plot      = p,
  save_path = save_name,
  width     = size$width,
  height    = size$height
)

# ───────────────────────────────────────────────────────────────
# Visualize: Correlation heatmap
# ───────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────
# Filter features based on technology (filtering done across all ancestries)
# ───────────────────────────────────────────────────────────────

filter_features  <- setup$filter_features
tech             <- setup$tech
percentile       <- setup$percentile
data_type        <- setup$data_type

save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")

if (filter_features && tech == "transcriptomics") {
  filtered_data <- filter_features_transcriptomics(
    adata      = adata,
    percentile = percentile,
    data_type  = data_type,
    plot_path  = save_name
  )
} else{
  stop("Filtering did not work!")
}

# ───────────────────────────────────────────────────────────────
# Ancestry Stratification
# ───────────────────────────────────────────────────────────────

save_name <- file.path(path_to_save_location, "QC_ancestry_stratification.pdf")

strat_result <- stratify_ancestry_subsets(
  adata           = filtered_data,
  ancestry_column = ancestry_column,
  output_column   = output_column,
  train_ancestry  = train_ancestry,
  infer_ancestry  = infer_ancestry,
  seed            = seed,
  plot_path       = save_name
)

train_adata <- strat_result$train_adata
test_adata  <- strat_result$test_adata
infer_adata <- strat_result$infer_adata

# ───────────────────────────────────────────────────────────────
# Calculate observed values (Non-Bootstrap Ground Truth)
# ───────────────────────────────────────────────────────────────

normalization_method <- setup$normalization_method
normalization_fn     <- match.fun(normalization_method)

train_obs <- calculate_logfc(
  matrix        = train_adata$X,
  meta          = train_adata$obs,
  group_column  = output_column,
  normalization = normalization_fn
)

test_obs <- calculate_logfc(
  matrix        = test_adata$X,
  meta          = test_adata$obs,
  group_column  = output_column,
  normalization = normalization_fn
)

infer_obs <- calculate_logfc(
  matrix        = infer_adata$X,
  meta          = infer_adata$obs,
  group_column  = output_column,
  normalization = normalization_fn
)

# ───────────────────────────────────────────────────────────────
# Calculate observed statistic
# ───────────────────────────────────────────────────────────────

interaction_obs <- compute_logfc_difference(
  test  = test_obs,
  infer = infer_obs
)

pearson_obs <- compute_logfc_correlation_difference(
  train  = train_obs,
  test   = test_obs,
  infer  = infer_obs,
  method = "pearson"
)

spearman_obs <- compute_logfc_correlation_difference(
  train  = train_obs,
  test   = test_obs,
  infer  = infer_obs,
  method = "spearman"
)

# ───────────────────────────────────────────────────────────────
# Bootstrapping
# ───────────────────────────────────────────────────────────────

H0_matrix <- rbind(test_adata$X, infer_adata$X)
H0_meta   <- rbind(test_adata$obs, infer_adata$obs)
n_boots   <- setup$n_bootstraps

test_boot <- bootstrap_logfc(
  matrix        = H0_matrix,
  meta          = H0_meta,
  size          = nrow(test_adata),
  group_column  = output_column,
  normalization = normalization_fn,
  n_iterations  = n_boots
)

infer_boot <- bootstrap_logfc(
  matrix        = H0_matrix,
  meta          = H0_meta,
  size          = nrow(infer_adata),
  group_column  = output_column,
  normalization = normalization_fn,
  n_iterations  = n_boots
)

# ───────────────────────────────────────────────────────────────
# Calculate bootstrapped statisitc
# ───────────────────────────────────────────────────────────────

interaction_boot <- compute_logfc_difference(
  test  = test_boot,
  infer = infer_boot
)

boot_ids <- unique(test_boot$bootstrap)
train_boot <- train_obs[
  , .(bootstrap = boot_ids), by = .(feature, logFC)
]

pearson_boot <- compute_logfc_correlation_difference(
  train  = train_boot,
  test   = test_boot,
  infer  = infer_boot,
  method = "pearson"
)

spearman_boot <- compute_logfc_correlation_difference(
  train  = train_boot,
  test   = test_boot,
  infer  = infer_boot,
  method = "spearman"
)

# ───────────────────────────────────────────────────────────────
# P-value
# ───────────────────────────────────────────────────────────────

save_name <- file.path(path_to_save_location, "QC_null_interaction.pdf")
interaction_summary <- calculate_pvalues(
  observed_dt  = interaction_obs,
  bootstrap_dt = interaction_boot,
  value_col    = "logFC",
  by_cols      = "feature",
  is_global    = FALSE,
  plot_path    = save_name,
  features     = interaction_obs[1:9, feature]
)

save_name <- file.path(path_to_save_location, "QC_null_pearson.pdf")
pearson_summary <- calculate_pvalues(
  observed_dt  = pearson_obs,
  bootstrap_dt = pearson_boot,
  value_col    = "difference",
  by_cols      = "feature",
  is_global    = TRUE,
  plot_path    = save_name,
  features     = NULL
)

save_name <- file.path(path_to_save_location, "QC_null_spearman.pdf")
spearman_summary <- calculate_pvalues(
  observed_dt  = spearman_obs,
  bootstrap_dt = spearman_boot,
  value_col    = "difference",
  by_cols      = "feature",
  is_global    = TRUE,
  plot_path    = save_name,
  features     = NULL
)


# ───────────────────────────────────────────────────────────────
# Adjused p-value
# ───────────────────────────────────────────────────────────────

interaction_summary[, p_adj_emp := p.adjust(p_emp, method = "BH")]
interaction_summary[, p_adj_param := p.adjust(p_param, method = "BH")]

# ───────────────────────────────────────────────────────────────
# Save the summary
# ───────────────────────────────────────────────────────────────

save_name <- file.path(path_to_save_location, "Interaction_summary.csv")
fwrite(interaction_summary, save_name)

correlation_summary <- rbind(
  x = pearson_summary,
  y = spearman_summary
)

save_name <- file.path(path_to_save_location, "Correlation_summary.csv")
fwrite(correlation_summary, save_name)

# ───────────────────────────────────────────────────────────────
# Visualization
# ───────────────────────────────────────────────────────────────

visual_val <- setup$visual_val

if (visual_val){

  path_to_save_location <- file.path(setup$output_directory, "Visual_val")
  if (!dir.exists(path_to_save_location)) {
    dir.create(
      path      = path_to_save_location,
      recursive = TRUE
    )
  }

  p <- plot_volcano(
    data         = interaction_summary,
    effect_col   = "logFC",
    p_value_col  = "p_adj_param",
    effect_label = "logFC",
    title        = NULL
  )

  save_name <- file.path(path_to_save_location, "Volcano_plot.pdf")

  size <- estimate_plot_size(p)

  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )

  p <- plot_correlation(
    data        = correlation_summary,
    test_label  = "Test",
    infer_label = "Inference"
    )
  
  save_name <- file.path(path_to_save_location, "Correlation_plot.pdf")

  size <- estimate_plot_size(p)

  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )


}

# ───────────────────────────────────────────────────────────────
# Sanity check: Limma interactions 
# ───────────────────────────────────────────────────────────────

sanity_check <- setup$sanity_check

if (sanity_check){

  path_to_save_location <- file.path(setup$output_directory, "Sanity_check_limma")
  if (!dir.exists(path_to_save_location)) {
    dir.create(
      path      = path_to_save_location,
      recursive = TRUE
    )
  }

  limma_matrix <- rbind(test_adata$X, infer_adata$X)
  limma_meta   <- rbind(test_adata$obs, infer_adata$obs)

  limma_res <- do_interactions(
    expr_data       = t(limma_matrix),
    sample_metadata = limma_meta,
    ancestry_column = ancestry_column,
    output_column   = output_column,
    train_ancestry  = train_ancestry,
    infer_ancestry  = infer_ancestry,
    comparison      = comparison
  )

  limma_means <- extract_top_table(
    limma_fit      = limma_res$fit_means,
    feature_column = "feature"
  )

  limma_contrasts <- extract_top_table(
    limma_fit      = limma_res$fit_contrasts,
    feature_column = "feature"
  )
  
  save_name <- file.path(path_to_save_location, "Limma_means.csv")
  fwrite(limma_means, save_name)

  save_name <- file.path(path_to_save_location, "Limma_contrast.csv")
  fwrite(limma_contrasts, save_name)

  limma_interactions <- limma_contrasts[grepl("Interaction", coef)]

  p <- plot_volcano(
    data         = limma_interactions,
    effect_col   = "logFC",
    p_value_col  = "adj.P.Val",
    effect_label = "logFC",
    title        = NULL
  )

  save_name <- file.path(path_to_save_location, "QC_volcano_plot.pdf")

  size <- estimate_plot_size(p)

  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )

}



