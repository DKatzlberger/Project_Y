#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────
# Cross-Ancestry Parametric Multi-Seed Script (Development Mode)
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

# ───────────────────────────────────────────────────────────────
# YAML settings file (from CLI or fallback)
# ───────────────────────────────────────────────────────────────

YAML_FILE <- yaml_argument_loader(
  default_yaml = "inst/extdata/example_settings_bootstrapped_interactions_parametric_multi_seed.yaml"
)

# ───────────────────────────────────────────────────────────────
# Default settings file
# ───────────────────────────────────────────────────────────────

DEFAULT_FILE <- get_default_settings_path(
  file_name = "default_settings_bootstrapped_interactions_parametric.yaml",
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
  "multi_seed",
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
# Set parallel execution strategy
# ───────────────────────────────────────────────────────────────

plan(multisession)

# ───────────────────────────────────────────────────────────────
# Run pipeline across seeds in parallel
# ───────────────────────────────────────────────────────────────

seeds <- setup$multi_seed
future_walk(
  seeds,
  function(seed_value) {
    cat(glue::glue("\n ▶ Running seed {seed_value}...\n"))

    # Load function files inside each worker
    invisible(lapply(function_files, source))

    args                  <- setup
    args$seed             <- seed_value
    args$output_directory <- file.path(setup$output_directory, paste0("seed_", seed_value))
    dir.create(args$output_directory, recursive = TRUE, showWarnings = FALSE)

    tryCatch(
      {
        do.call(bootstrapped_interactions_parametric_pipeline, args)
      },
      error = function(e) {
        error_log <- file.path(args$output_directory, "error_log.txt")
        msg <- glue::glue("
          SEED {seed_value} FAILED
          ────────────────────────────────────────────
          Message: {conditionMessage(e)}
          Call:    {deparse(conditionCall(e))}
          ────────────────────────────────────────────
          Traceback:
          {paste(capture.output(traceback(max.lines = 5)), collapse = '\n')}
          ")
        cat(msg, file = stderr())             # Print to console
        writeLines(msg, con = error_log)      # Save to file
      }
    )
  },
  .options = furrr_options(seed = TRUE)
)
