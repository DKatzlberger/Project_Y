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
analysis_name = "visualization"

tag <- setup$tag
match_pattern <- paste0(tag, "_", analysis_name)
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Load the data
adata <- read_h5ad(setup$data_path)


# Descriptive statistics
remove <- c("sample_id", 
            "study_id", 
            "name", 
            "cancer_type_id", 
            "cancer_type_name",
            "patient_id",
            "tumor_type",
            "rnaseqv2"
            )

remove <- c(remove)
# Remove columns with unnessary information           
to_visualize <- setdiff(colnames(adata$obs), remove)
# Classify meta data
calssified_meta <- classify_meta(adata$obs, to_visualize)

# Visualize: Discrete variables
discrete_plots <- lapply(calssified_meta$discrete, function(var) {
  ggplot(adata$obs, aes(x = !!sym(var))) +
    geom_bar() +
    facet_grid(rows = vars(!!sym(setup$classification$ancestry_column)), scales = "free") +
    labs(title = var, x = var, y = "Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5))
})
# Combine
all_discrete <- wrap_plots(discrete_plots, ncol = max(floor(length(discrete_plots) / 2), 1))
# Calculate dimensions for figure
n_cols = max(floor(length(discrete_plots) / 2), 1)
n_rows = length(discrete_plots) / n_cols
width = n_cols * 4 
height = n_rows * 7
# Save
ggsave(filename = "Discrete_variables.pdf", 
       plot = all_discrete, 
       path = path_to_save_location, 
       width = width, height = height
       )

# Visualize: Continous variables
continuous_plots <- lapply(calssified_meta$continuous, function(var) {
  ggplot(adata$obs, aes(x = !!sym(var))) +
    geom_density() +
    facet_grid(rows = vars(!!sym(setup$classification$ancestry_column)), scales = "free") +
    labs(title = var, x = var, y = "Density") +
    theme(plot.title = element_text(hjust = 0.5))
})
# Combine
all_continous <- wrap_plots(continuous_plots, ncol = max(floor(length(continuous_plots) / 2), 1))
# Calculate dimensions for figure
n_cols = max(floor(length(continuous_plots) / 2), 1)
n_rows = length(continuous_plots) / n_cols
width = n_cols * 4
height = n_rows * 7
# Save
ggsave(filename = "Continous_variables.pdf", 
       plot = all_continous, 
       path = path_to_save_location, 
       width = width, height = height
       )

# Clustering analysis
# Check NA values
# Rows
rows_with_any_na_count <- sum(apply(adata$X, 1, function(row) any(is.na(row))))
rows_with_all_na_count <- sum(apply(adata$X, 1, function(row) all(is.na(row))))
# Cols
columns_with_any_na_count <- sum(apply(adata$X, 2, function(col) any(is.na(col))))
columns_with_all_na_count <- sum(apply(adata$X, 2, function(col) all(is.na(col))))

# Remove features (columns) that contain NA
counts_na_removed <- adata$X[, apply(adata$X, 2, function(col) !any(is.na(col)))]
# Check variance
counts_with_variance <- counts_na_removed[, apply(counts_na_removed, 2, var) > 0]
# Transform  data to log2-normalized counts
log2_normalized <- log2(counts_with_variance + 1)

# PCA
pca_result <- prcomp(log2_normalized, scale. = TRUE)
# Filter first 5 PCs
filtered_pcs <- pca_result$x |>
  as_tibble(rownames = "patient_id") |>
  select(patient_id, PC1, PC2, PC3, PC4, PC5, PC6) |>
  merge(as_tibble(adata$obs, rownames = "patient_id"), by = "patient_id")

# Plot colored by each discrete variable
for (condition in calssified_meta$discrete){
  # Generate PCA plots for the condition
  pca_plots <- plot_pca_for_condition(filtered_pcs, 
                                      pca_result = pca_result, 
                                      condition = condition
                                      )
  # Combine the plots with patchwork
  combined_plot <- wrap_plots(pca_plots, ncol = 4) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",           
        legend.direction = "horizontal"      
      )
  # Save the combined plot
  ggsave(
    filename = paste0("PCA_", condition, ".pdf"), 
    plot = combined_plot, 
    path = path_to_save_location, 
    width = 11, height = 11
  )
}

# TSNE
# Reduce to 2 dimensions
perplexities = round(seq(5, 50, length.out = 9))
for (condition in calssified_meta$discrete){
  # Generate T-SNE plots for condition
  tsne_plots <- generate_tsne_plots(data = log2_normalized,
                                    meta = adata$obs,
                                    condition = condition,
                                    perplexity_range = perplexities
                                    )
  # Combine plots
  combined_plot <- wrap_plots(tsne_plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",           
      legend.direction = "horizontal"      
      )
  # Save 
  ggsave(
    filename = paste0("T-SNE_", condition, ".pdf"), 
    plot = combined_plot, 
    path = path_to_save_location, 
    width = 11, height = 11
  )
}