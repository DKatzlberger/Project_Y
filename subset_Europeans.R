# Remove start up messages
suppressPackageStartupMessages(
  {
    # Standard libraries
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # Statistics library
    library(edgeR)
    # Visualization
    library(patchwork)
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
  is_yml_file(YAML_FILE)
  setup <- yaml.load_file(YAML_FILE)

} else {
  # Dev settings if no command-line argument provided
  yaml_file <- "job_settings.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}

# Create path to output directory
# base_dir = "data/combined_runs"
# tag <- setup$tag
# comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
# analysis_name = "Interactions"

# # Create directory if not exist
# output_path = file.path(base_dir, comparison, analysis_name)
# if (!dir.exists(output_path)) {
#   dir.create(output_path, recursive = TRUE)
# }

# Test on Europeans how correlation changes with decrease in sample size
# This is done for each classification task
# Split all Europeans 50-50 but randomly, this means you will use seeds
# I do not stratify the samples

# Load the data
adata <- read_h5ad(setup$data_path)
# Set seed (Split of EUR is done randomly)
seed <- setup$seed
set.seed(setup$seed)

# Define classification task
data <- adata[adata$obs[[setup$classification$output_column]]
              %in% setup$classification$comparison]

# Filter only EUR based on ancestry column (column of used ancestry)
eur_data <- data[data$obs[[setup$classification$ancestry_column]] 
                 == setup$classification$train_ancestry]

# Randomly split EUR into two 50 50 subsets (currently I dont stratify)
# Generate a random index to split the data frame
split_index <- sample(1:nrow(eur_data), nrow(eur_data) / 2)
# Create two EUR subsets based on the random index
stable_subset <- eur_data[split_index, ]
stable_n_obs <- nrow(stable_subset)

sampling_subset <- eur_data[-split_index, ]
# The 'stable_subset' is keept as is (no further splitting)
# The 'sampling_subset' is used for downsampling and stored in a list of subsets

# For the 'stable_subset' logFC are obtained using limma
# Creating design matrix for the 
stable_design <- create_design(setup$classification$output_column,
                               meta = stable_subset$obs)

# Create contrast matrix (only needs to be created once)
# Used for hypothesis testing between groups
contrast_matrix <- create_contrast(colnames(stable_design))

# Filter lowley expressed genes bases on 'stable_subset'
keeper_genes <- filterByExpr(t(stable_subset$X), stable_design)
# Remove lowley expressed genes
stable_subset <- stable_subset[, keeper_genes]
sampling_subset <- sampling_subset[, keeper_genes]

# Limma for the 'stable_subset'
# logCPM transformation
stable_norm <- voom(t(stable_subset$X), stable_design, plot = FALSE)

# Fit the model (means model)
stable_limma_fit <- lmFit(stable_norm, design = stable_design)
# Ebayes
stable_limma_fit <- eBayes(stable_limma_fit)
# Results
stable_means_res <- extract_results(stable_limma_fit)

# Fit contrast
stable_limma_fit_contrast <- contrasts.fit(stable_limma_fit, contrast_matrix)
# Ebayes
stable_limma_fit_contrast <- eBayes(stable_limma_fit_contrast)
# Results
stable_contrast_res <- extract_results(stable_limma_fit_contrast)

# Get logFC for correlation analysis 
stable_logFC <- stable_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_stable = logFC) 

# After filtering of genes take samples from 'sampling_subset'
# This is not stratified!
props <- c(1.0, 0.8, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)
# In order for eBayes replicates are needed 
# Check that samples have replicates per cancer
sample_list <- sample_with_proportions(sampling_subset, 
                                       proportions = props, 
                                       seed = setup$seed)

# Visualization of the distribution of obervation in the subsets
# Currently the distribution of the classes is random
sample_meta <- data.frame()
for (i in seq_along(sample_list)) {
  smpl <- sample_list[[i]]
  sample_name <- names(sample_list)[i]
  # Get meta and add column Sample
  meta <- smpl$obs |> 
    mutate(Sample = sample_name)
  # Combine
  sample_meta <- bind_rows(sample_meta, meta) 
}
# Add meta of 'stable_subset'
sample_meta <- 
  stable_subset$obs|>
  mutate(Sample = "stable_subset") |>
  bind_rows(sample_meta)

# Plot with ggplot
sample_distribution <- ggplot(
    sample_meta,
    aes(
        x = fct_rev(Sample),
        fill = get(setup$classification$output_column)
        )
    ) +
    geom_bar(
        position = "stack"
    ) +
    theme(axis.text.x = element_text(angle = 90))


# Limma for all other subsets
logFC_results <- NULL
# Iterate over all samples in sample_list
for (i in seq_along(sample_list)) {

  subset <- sample_list[[i]]
  # Check if observations in subset are unique
  stopifnot(anyDuplicated(subset$obs_names) == 0)

  # Information about the subset
  subset_name <- names(sample_list)[i]
  n_obs <- nrow(subset)
  # Create design matrix
  subset_design <- create_design(setup$classification$output_column, meta = subset$obs)

  # logCPM transformation
  subset_norm <- voom(t(subset$X), subset_design, plot = FALSE)

  # Fit the model (means model)
  subset_limma_fit <- lmFit(subset_norm, design = subset_design)
  # Ebayes
  subset_limma_fit <- eBayes(subset_limma_fit)
  # Results
  susbet_means_res <- extract_results(subset_limma_fit)

  # Fit contrast
  subset_limma_fit_contrast <- contrasts.fit(subset_limma_fit, contrast_matrix)
  # Ebayes
  subset_limma_fit_contrast <- eBayes(subset_limma_fit_contrast)
  # Results
  subset_contrast_res <- extract_results(subset_limma_fit_contrast)

  # LogFC for correlation analysis
  subset_logFC <- subset_contrast_res |>
    select(coef, Feature, logFC) |>
    rename_with(~ paste0("logFC_", subset_name), .cols = logFC) 
  
  # Combine results
  if (is.null(logFC_results)) {
    # Initialize with the first subset
    logFC_results <- subset_logFC
  } else {
    # Merge with existing results
    logFC_results <- inner_join(logFC_results, subset_logFC, by = c("coef", "Feature"))
  }
}

# Add the logFC of the 'stable_subset' (ground truth)
logFC_results <- inner_join(logFC_results, stable_logFC, by = c("coef", "Feature"))

# Correlation of all numeric columns (logFC)
pearson <- compute_correlation(data = logFC_results, method = "pearson") |>
  filter(str_detect(V2, "stable"))
spearman <- compute_correlation(data = logFC_results, method = "spearman") |>
  filter(str_detect(V2, "stable"))

# Make metric 

