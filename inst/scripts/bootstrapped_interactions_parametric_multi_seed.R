# Script to run cross-ancestry pipeline across multiple seeds in parallel

library(furrr)
library(purrr)

# Load your package or source the pipeline function if not yet in a package
# library(yourPackageName)
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

# Set parallel execution strategy
plan(multisession)

# Define seeds
seeds <- c(1, 42, 43, 44)

# Define base arguments (without seed and output_directory)
base_args <- list(
  output_column        = "subtype_pooled",
  class_0              = "Basal",
  class_1              = "non-Basal",
  ancestry_column      = "genetic_ancestry",
  train_ancestry       = "eur",
  infer_ancestry       = "afr",
  data_path            = "data/inputs/PanCanAtlas_BRCA_raw_RSEM_subtypeNAremoved.h5ad",
  tech                 = "transcriptomics",
  normalization_method = "voom_normalization",
  filter_features      = TRUE,
  percentile           = 25,
  data_type            = "RSEM",
  n_bootstraps         = 1000,
  sanity_check         = TRUE,
  save_outputs         = TRUE
)

# Run the pipeline for each seed in parallel with safe random seeding
future_walk(
  seeds,
  function(s) {
    output_dir <- file.path("results/multi_seed", paste0("seed_", s))

    args <- base_args
    args$seed <- s
    args$output_directory <- output_dir

    if (is.null(output_dir) || output_dir == "") {
      stop("Output directory is not defined.")
    }

    do.call(bootstrapped_interactions_parametric_pipeline, args)
  },
  .options = furrr_options(seed = TRUE)
)
