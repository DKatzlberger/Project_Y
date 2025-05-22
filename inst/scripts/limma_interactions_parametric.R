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
cat("-----------------------------------------------------------\n")

# ───────────────────────────────────────────────────────────────
# YAML settings file (from CLI or fallback)
# ───────────────────────────────────────────────────────────────

YAML_FILE <- yaml_argument_loader(
  default_yaml = "inst/extdata/example_settings_interactions_parametric.yaml"
)

# ───────────────────────────────────────────────────────────────
# Merge default + user YAML settings
# ───────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────
# Default settings file
# ───────────────────────────────────────────────────────────────

DEFAULT_FILE <- get_default_settings_path(
  file_name = "default_settings_interactions_parametric.yaml",
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
# Interactions
# ───────────────────────────────────────────────────────────────

limma_res <- do_interactions(
  expr_data       = t(filtered_data$X),
  sample_metadata = filtered_data$obs,
  ancestry_column = ancestry_column,
  output_column   = output_column,
  train_ancestry  = train_ancestry,
  infer_ancestry  = infer_ancestry,
  comparison      = comparison
)

means <- extract_top_table(
  limma_fit      = limma_res$fit_means,
  feature_column = "feature"
)

contrasts <- extract_top_table(
  limma_fit      = limma_res$fit_contrasts,
  feature_column = "feature"
  )

# ───────────────────────────────────────────────────────────────
# Save the results
# ───────────────────────────────────────────────────────────────

save_name <- file.path(path_to_save_location, "Limma_means.csv")
fwrite(means, save_name)

save_name <- file.path(path_to_save_location, "Limma_contrast.csv")
fwrite(contrasts, save_name)