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
print("Evoking R...")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  yaml_file <- args[1]
  # Check if it's a valid YAML file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
} else {
  print("Running interactive mode for development.")
  # Dev settings if no command-line argument provided
  yaml_file <- "dev_settings.yml"
  setup <- yaml.load_file(yaml_file)
}

# Set seed (because why not)
set.seed(setup$seed)

# Load data
adata <- read_h5ad(setup$data_path)

# Load the selected feature used in ml
output_directory <- setup$output_directory
features <- yaml.load_file(file.path(output_directory, "Features.yml"))
feature_n <- length(features)

# Start with 'constant_split'
print("Dealing with constant split.")
# Load and subset
constant_idx <- yaml.load_file(file.path(output_directory, "Obs_constant.yml"))
constant_split <- adata[constant_idx, features]
constant_n = length(constant_idx)

# Assertion: Check arrays
# Dimensions
stopifnot(length(constant_idx) == length(constant_split$obs_names))
# Order
stopifnot(all(colnames(t(constant_split$X)) == row.names(constant_split$obs[colnames(t(constant_split$X)), ])))

# Covariate (needs to be treated differently)
# 1. Extract the covariate, if it exists and has a meaningful value
covariate_list <- get_covariate(setup)
# 2. Check if there are replicates and at least two levels for the covariate 
if (!is.null(covariate_list)){
  # Check covariate type
  print("Checking covariates.")
  covariate_types <- classify_covariates(constant_split$obs, covariate_list)
  # Check conditions "Discrete"
  check_covariate_conditions(constant_split$obs, covariates = covariate_types$discrete)
}

# Create design matrix
output_column <- setup$classification$output_column
constant_design <- create_design(output_column, constant_split$obs, covariate = covariate_list)

# Limma workflow
print("Start differential gene expression analysis.")
# Select normalization method
data_type <- setup$data_type
dge_normalization <- setup$dge_normalization
normalization_method <- normalization_methods[[data_type]][[dge_normalization]]

# Transpose (rows = Genes, cols = Samples)
constant_split_t <- t(constant_split$X)

# Normalization
constant_norm <- normalization_method(constant_split_t, constant_design)

# Fit the model (means model)
constant_limma_fit <- lmFit(constant_norm, design = constant_design)
constant_limma_fit <- eBayes(constant_limma_fit)
constant_means_res <- extract_results(constant_limma_fit)

# Create contrast matrix 
print("Fit contrast.")
comparison  <- setup$classification$comparison
contrast_matrix_constant <- create_contrast(colnames(constant_design), conditions = comparison)
# Fit contrast
constant_limma_fit_contrast <- contrasts.fit(constant_limma_fit, contrast_matrix_constant)
constant_limma_fit_contrast <- eBayes(constant_limma_fit_contrast)
constant_contrast_res <- extract_results(constant_limma_fit_contrast)
# Results constant split
constant_logFC <- constant_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_constant = logFC) 

# Subsets
vscratch_dir <- setup$output_directory
match_pattern <- paste0("Obs_proportion_")

# Extracting all folders in the 'vscratch_dir' that match 'match_pattern'
all_vscratch_dir <- list.files(vscratch_dir, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir, value = TRUE)

# Create a list with all subsets
subset_list <- list()
for (file in match_vscratch_dir){
    # Filter proportion
    proportion_pattern <- gsub("Obs_proportion_", "", basename(file))  
    proportion_pattern <- gsub(".yml", "", proportion_pattern)        
    proportion_pattern <- gsub("_", ".", proportion_pattern)  
    # Observations
    index <- yaml.load_file(file)
    subset <- adata[index, features]
    # Append to list
   subset_list[[proportion_pattern]] <- subset
}

# Limma analysis for all other subsets
logFC_results <- NULL
for (i in seq_along(subset_list)){

  prop <- names(subset_list)[[i]]
  print(paste("Dealing with", prop, "of data."))

  subset <- subset_list[[i]]
  # Check if observations in subset are unique
  stopifnot(anyDuplicated(subset$obs_names) == 0)
  # Information about the subset
  n_obs <- nrow(subset)

  # Covariate
  if (!is.null(covariate_list)){
    # Check covariate type
    print("Checking covariates.")
    covariate_types <- classify_covariates(subset$obs, covariate_list)
    # Check conditions "Discrete"
    check_covariate_conditions(subset$obs, covariates = covariate_types$discrete)
  }

  # Design
  output_column <- setup$classification$output_column
  subset_design <- create_design(output_column, meta = subset$obs, covariate = covariate_list)
  # Transpose
  subset_t <- t(subset$X)
  # Normalization
  subset_norm <- normalization_method(subset_t, design = subset_design)
  # Limma
  subset_limma_fit <- lmFit(subset_norm, design = subset_design)
  subset_limma_fit <- eBayes(subset_limma_fit)
  susbet_means_res <- extract_results(subset_limma_fit)
  # Contrast
  contrast_matrix_subset <- create_contrast(colnames(subset_design), conditions = comparison)
  # Fit contrast
  subset_limma_fit_contrast <- contrasts.fit(subset_limma_fit, contrast_matrix_subset)
  subset_limma_fit_contrast <- eBayes(subset_limma_fit_contrast)
  subset_contrast_res <- extract_results(subset_limma_fit_contrast)
  # Results
  subset_logFC <- subset_contrast_res |>
    select(coef, Feature, logFC) |>
    rename_with(~ paste0("logFC_", n_obs), .cols = logFC) 

  # Combine results
  if (is.null(logFC_results)) {
    # Initialize with the first subset
    logFC_results <- subset_logFC
  } else {
    # Merge with existing results
    logFC_results <- inner_join(logFC_results, subset_logFC, by = c("coef", "Feature"))
  }

}

# Add the logFC of the 'constant_split' (ground truth)
logFC_results <- inner_join(logFC_results, constant_logFC, by = c("coef", "Feature"))

# Per gene metric
raw_logFC <- logFC_results |>
  pivot_longer(
    cols = starts_with("logFC"),
    names_to = "n_test_ancestry", 
    values_to = "logFC"
  ) |>
  mutate(
    n_test_ancestry = str_remove(n_test_ancestry, "logFC_")
  ) |>
  mutate(
      Seed = setup$seed,
      Ancestry = toupper(setup$classification$train_ancestry),
      n_train_ancestry = constant_n,
      Status = ifelse(n_test_ancestry == "constant", "Train", "Test")
  )

# Save
fwrite(raw_logFC, file.path(setup$output_directory, "LogFCs.csv"))

# Summarized metric
pearson <- compute_correlation(data = logFC_results, method = "pearson") |>
  filter(str_detect(V2, "constant"))
spearman <- compute_correlation(data = logFC_results, method = "spearman") |>
  filter(str_detect(V2, "constant"))

# Combine Pearson and Spearman
# Add all information need to combine multiple seeds
metric_df <- inner_join(pearson, spearman, by = c("V1", "V2")) |>
  mutate(
    Ancestry = toupper(setup$classification$train_ancestry),
    Seed = setup$seed,
    n_test_ancestry = str_extract(V1, "\\d+$"),
    n_train_ancestry = constant_n
  )

# Save metric Dataframe
fwrite(metric_df, file.path(setup$output_directory, "Contrast_metric_dge.csv"))


if (setup$visual_val){
  # Plot all subsets
  # 1. Iterate over all subsets
  # 2. Collect meta for sampled subsets
  # 3. Add 'constant_split' meta
  subset_meta_combined <- data.frame()
  for (i in seq_along(subset_list)){
    subset <- subset_list[[i]]
    subset_meta <- subset$obs
    # Add additional information
    subset_meta <- subset_meta |> mutate(n_obs = n(), 
                                        type = "sampled"
                                        )
    # Append
    subset_meta_combined <- bind_rows(subset_meta_combined, subset_meta)
  }

  # 'constant_split'
  constant_meta <- constant_split$obs
  # Add additional information
  constant_meta <- constant_meta |> mutate(n_obs = n(), 
                                          type = "constant"
                                          )
  # Append to 'subset_meta_combined'
  subset_meta_combined <- bind_rows(subset_meta_combined, constant_meta)

  # Visualite meta data
  # Bar plot
  subset_plot <- 
    ggplot(
      subset_meta_combined,
      aes(
        x = fct_rev(as_factor(n_obs)),
        fill = get(setup$classification$output_column)
      )
    ) +
    geom_bar(
      position = "stack"
    ) +
    facet_grid(
      cols = vars(type),
      scales = "free",
      space = "free"
    ) +
    labs(
      x = "Number of observation",
      y = "Count",
      fill = setup$classification$output_column,
    ) + 
    theme(
      plot.caption = element_text(hjust = 0)
    ) 

  # Save 'subset_plot'
  ggsave(file.path(setup$output_directory, "Validation_stratification.pdf"),
        plot = subset_plot, height = 6, width = 12)
}

# Print statement to switch back to python script
print("Switching back to Python.")
