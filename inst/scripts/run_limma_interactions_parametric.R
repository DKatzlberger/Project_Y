#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────
# Interactions Parametric Main Script (Development Mode)
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(anndata)
  library(yaml)
  library(glue)
  library(uuid)
  library(reticulate)
  library(parallel)
  library(furrr)
  library(limma)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)

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
  default_yaml = "inst/extdata/example_settings_limma_interactions_parametric.yaml"
)

# ───────────────────────────────────────────────────────────────
# Default settings file
# ───────────────────────────────────────────────────────────────

DEFAULT_FILE <- get_default_settings_path(
  file_name = "default_settings_limma_interactions_parametric.yaml",
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
  # Classification
  "output_column", 
  "class_0", 
  "class_1", 
  "ancestry_column", 
  "train_ancestry", 
  "infer_ancestry", 
  "data_path", 
  "tech", 
  "output_directory",
  "normalization_method"
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

save_name <- file.path(path_to_save_location, "QC_correaltion_heatmap.pdf")

plot_correlation_heatmap(
  expr_data       = t(adata$X),
  metadata        = adata$obs,
  ancestry_column = ancestry_column,
  output_column   = output_column,
  file_path       = save_name,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_names      = FALSE,
  width           = 5,
  height          = 5
)


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

# Save features 
save_name <- file.path(path_to_save_location, "Features.yaml")
write_yaml(filtered_data$var_names, save_name)

setup$n_features <- length(filtered_data$var_names)
save_settings(
  settings   = setup,
  file_path  = file.path(setup$output_directory, "Settings.yaml")
)

# ───────────────────────────────────────────────────────────────
# Interactions
# ───────────────────────────────────────────────────────────────

limma_res <- do_limma_interactions(
  expr_data       = t(filtered_data$X),
  sample_metadata = filtered_data$obs,
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

# ───────────────────────────────────────────────────────────────
# Save the results
# ───────────────────────────────────────────────────────────────

save_name <- file.path(path_to_save_location, "Limma_means.csv")
fwrite(limma_means, save_name)

save_name <- file.path(path_to_save_location, "Limma_contrast.csv")
fwrite(limma_contrasts, save_name)

limma_interactions <- limma_contrasts[grepl("Interaction", coef)]
limma_sig_features <- limma_interactions[abs(logFC) > 1 & adj.P.Val < 0.05]$feature

save_name <- file.path(path_to_save_location, "Significant_features.yaml")
write_yaml(limma_sig_features, save_name)

# ───────────────────────────────────────────────────────────────
# Visualize volcano
# ───────────────────────────────────────────────────────────────

path_to_save_location <- file.path(path_to_save_location, "Visual_val")
if (!dir.exists(path_to_save_location)) {
  dir.create(
    path      = path_to_save_location,
    recursive = TRUE
  )
}

p <- plot_volcano(
  data         = limma_interactions,
  effect_col   = "logFC",
  p_value_col  = "adj.P.Val",
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

# ───────────────────────────────────────────────────────────────
# Visualize significant genes (Heatmap/Boxplots)
# ───────────────────────────────────────────────────────────────

if (lenght(limma_sig_features) != 0){

  zscores <- prepare_normalized_data(
    matrix   = limma_res$norm_matrix,
    meta     = adata$obs,
    features = limma_sig_features
    )
  
  # Boxplots
  p <- plot_expression_boxplots(
    data = zscores,
    x    = ancestry_column,
    y    = "zscore",
    fill = output_column
  )

  save_name <- file.path(path_to_save_location, "Interaction_boxplot.pdf")

  size <- estimate_plot_size(p)

  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )

  # Expression heatmap
  save_name <- file.path(path_to_save_location, "Interaction_heatmap.pdf")
  p <- plot_expression_heatmap(
    zscores       = zscores,
    ancestry_column = ancestry_column,
    output_column   = output_column,
    file_path       = save_name,
    width           = 10,
    height          = 5
  )

}