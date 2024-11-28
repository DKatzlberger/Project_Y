# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Statistics 
    library(coin)
    }
)
# Source custom functions
source("r_utils.R")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  yaml_file <- args[1]  # First argument is the YAML file path
} else {
  # Dev settings
  yaml_file <- "job_settings.yml"
  print("Running interactive mode for development.\n")
}

# Load setup file
setup <- yaml.load_file(yaml_file)
# Location where I store my files
vscratch_dir <- "data/combined_runs"


# Tag is used to know which data is used
tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
train_ancestry <- toupper(setup$classification$train_ancestry)
infer_ancestry <- toupper(setup$classification$infer_ancestry)

# Directory for 
directory <- file.path(vscratch_dir, comparison)
all_folders <- list.dirs(directory, full.names = TRUE, recursive = FALSE)
matching_folders <- grep(comparison, all_folders, value = TRUE)

for (folder in matching_folders){
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}


