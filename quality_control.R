# Remove start up messages
suppressPackageStartupMessages(
  {
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # DGE workflow and functional analysis
    library(Rtsne)

    # Visualization
    library(patchwork)
    library(ggrepel)
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
  }
)
# Custom functions
source("r_utils.R")

# Here starts the script

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  yaml_file <- args[1]
  # Check if it's a valid YAML file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

} else {
  # Dev settings if no command-line argument provided
  yaml_file <- "dev_settings.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}
# Set seed (needed for tSNE)
set.seed(42)

# Construction:
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out = "data/combined_runs"
analysis_name = "EUR_to_EAS_quality_control"

tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
match_pattern <- paste0(comparison, "_", analysis_name)
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Load the data
adata <- read_h5ad(setup$data_path)

# Define classification task
data <- adata[adata$obs[[setup$classification$output_column]]
              %in% setup$classification$comparison]

# Filter by ancestry 
to_analyse_ancestries <- c(setup$classification$train_ancestry,
                           setup$classification$infer_ancestry)

data <- data[data$obs[[setup$classification$ancestry_column]]
             %in% to_analyse_ancestries]

# Visualize meta
remove <- c("sample_id", 
            "study_id", 
            "name", 
            "cancer_type_id", 
            "cancer_type_name",
            "patient_id"
            )
# Remove columns with unnessary information           
to_visualize <- setdiff(colnames(data$obs), remove)
# Classify meta data
calssified_meta <- classify_meta(data$obs, to_visualize)

discrete_plots <- lapply(calssified_meta$discrete, function(var) {
  ggplot(data, aes_string(x = var)) +
    geom_bar() +
    labs(title = paste("Bar plot of", var)) +
    theme_minimal()
})
