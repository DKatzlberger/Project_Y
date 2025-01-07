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
  yaml_file <- args[1] 
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
  # When the script is run from the command line then 'output_directory' is given
  # The pattern to extract all matchin directories is extracted from 'output_directory'
  # output_path = setup$output_directory
  # match_pattern <- gsub("_\\d+$", "", output_path)
  # TODO - Check if regular expression also works for other approaches

} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "job_settings.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
}

# Construction:
# Vscratch_dir is the place where the files are stored
vscratch_dir_in = file.path("data", "runs")
# Tag is used to specify which data it is e.g. TCGA, NIAGADS
tag <- setup$tag
# Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
# Analysis_name specifies what was analysed
# E.g. comparsion of ancestries EUR_to_AMR, subsetting_EUR, robustness_AFR
# This often is modified depending which analysis
infer_ancestry <- toupper(setup$classification$infer_ancestry)
analysis_name  <- paste0(train_ancestry, "_to_", infer_ancestry, "_robustness")
# Combination of components to create the 'match_pattern'
# The 'match_pattern' is used as pattern to extract all folders in the vscratch dir
match_pattern <- paste0(comparison, "_", analysis_name)

# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vsratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Check if there were matching folders
if (length(match_vsratch_dir) == 0) {
  message("No matching folders found.")
}

# Save the results of the analysis
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, comparison, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}