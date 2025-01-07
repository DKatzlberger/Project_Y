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
print("Evoking R.")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if script is run with a command-line argument
if (length(args) > 0) {
  yaml_file <- args[1]
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "dev_settings.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
}

print("Start differential gene expression analysis.")

# Set seed (because why not)
set.seed(setup$seed)

# Load data
adata <- read_h5ad(setup$data_path)

# Get necessary information to load the subsest created
# Proportions
proportions = setup$proportion
# Save directory
vscratch_dir <- setup$output_directory

# Iterate over all proportions doing DGE analysis
robustness_results <- list()
for (prop in proportions){
  # Replace . with _ in 'prop'
  prop <- gsub("\\.", "_", prop)

  # Create file names
  observation_file_train <- paste0("Obs_proportion_", prop, "_train.yml")
  observation_file_test <- paste0("Obs_proportion_", prop, "_test.yml")
  observation_file_inf <- paste0("Obs_proportion_", prop, "_inf.yml")

  # Load indeces
  train_idx <- yaml.load_file(file.path(setup$output_directory, observation_file_train))
  test_idx <- yaml.load_file(file.path(setup$output_directory, observation_file_test))
  inf_idx <- yaml.load_file(file.path(setup$output_directory, observation_file_inf))

  # Get number of samples
  train_n = length(train_idx)
  test_n = length(test_idx)
  inf_n = length(inf_idx)

  # Subset the data by the indexes create in python
  train_data <- adata[adata$obs_names %in% train_idx, ]
  test_data <- adata[adata$obs_names %in% test_idx, ]
  inf_data <- adata[adata$obs_names %in% inf_idx, ]

  # Assertion: Check arrays
  # Dimensions
  stopifnot(length(train_idx) == length(train_data$obs_names))
  # Order
  stopifnot(all(colnames(t(train_data$X)) ==
                  row.names(train_data$obs[colnames(t(train_data$X)), ])))


  # Covariate
  # 1. Extract the covariate, if it exists and has a meaningful value
  # 2. Check if covariate is discrete or continues
  covariate <- if ("covariate" %in% names(setup$classification) && 
                  !is.null(setup$classification$covariate) &&
                  setup$classification$covariate != "") {
    setup$classification$covariate
  } else {
    NULL
  }

  if (!is.null(covariate)){
    # TODO - When age is available check if is continous
    covariate_type <- check_covariate_type(train_data, covariate = covariate)

    # Check conditions for 'covariate_type'
    # Discrete:
    if (covariate_type == "Discrete"){
      # Possible values of the covariate
      possible_values <- unique(train_data$obs[[covariate]])
      # Check conditions
      # 1. Check values of covariate
      # 2. Ensure that there are 2 samples per value
      # 3. The script gets terminated if conditions are not met
      check_covariate_conditions(train_data, covariate = covariate, possible_values = possible_values)
      check_covariate_conditions(test_data, covariate = covariate, possible_values = possible_values)
      
      # 'inf_data' doesn't change per seed
      # If conditions don't met they will never met
      check_covariate_conditions(inf_data, covariate = covariate, possible_values = possible_values)
    }
  }

  # Create design matrix
  train_design <- create_design(setup$classification$output_column,
                                meta = train_data$obs,
                                covariate = covariate)                            
  test_design <- create_design(setup$classification$output_column,
                              meta = test_data$obs,
                              covariate = covariate)
  inf_design <- create_design(setup$classification$output_column,
                              meta = inf_data$obs,
                              covariate = covariate)

  # Create contrast matrix (Only needs to be created once)
  # Used for hypothesis testing between groups
  contrast_matrix <- create_contrast(colnames(train_design))

  # Filter genes by expression (Use the train set)
  # This will lead to different features across proportions
  # TODO - This will lead to different features across proportions
  keeper_genes <- filterByExpr(t(train_data$X), train_design)
  # Subset
  train_filtered <- train_data[, keeper_genes]
  test_filtered <- test_data[, keeper_genes]
  inf_filtered <- inf_data[, keeper_genes]

  # Save features (for use in ML)
  # Features are saved per 'prop'
  feature_file <- paste0("Features_", prop, ".yml")
  write_yaml(train_filtered$var_names, file.path(setup$output_directory, feature_file))

  # Limma workflow
  # Normalization (logCPM)
  train_norm <- voom(t(train_filtered$X), train_design, plot = FALSE)
  test_norm <- voom(t(test_filtered$X), test_design, plot = FALSE)
  inf_norm <- voom(t(inf_filtered$X), inf_design, plot = FALSE)

  # Fit the model (Means model)
  train_limma_fit <- lmFit(train_norm, design = train_design)
  test_limma_fit <- lmFit(test_norm, design = test_design)
  inf_limma_fit <- lmFit(inf_norm, design = inf_design)

  # Ebayes
  train_limma_fit <- eBayes(train_limma_fit)
  test_limma_fit <- eBayes(test_limma_fit)
  inf_limma_fit <- eBayes(inf_limma_fit)

  # Results means model
  train_mean_res <- extract_results(train_limma_fit)
  test_mean_res <- extract_results(test_limma_fit)
  inf_mean_res <- extract_results(inf_limma_fit)

  # Fit contrast
  train_limma_fit_contrast <- contrasts.fit(train_limma_fit, contrast_matrix)
  test_limma_fit_contrast <- contrasts.fit(test_limma_fit, contrast_matrix)
  inf_limma_fit_contrast <- contrasts.fit(inf_limma_fit, contrast_matrix)

  # Ebayes
  train_limma_fit_contrast <- eBayes(train_limma_fit_contrast)
  test_limma_fit_contrast <- eBayes(test_limma_fit_contrast)
  inf_limma_fit_contrast <- eBayes(inf_limma_fit_contrast)

  # Results contrast
  train_contrast_res <- extract_results(train_limma_fit_contrast)
  test_contrast_res <- extract_results(test_limma_fit_contrast)
  inf_contrast_res <- extract_results(inf_limma_fit_contrast)

  # Save contrast results
  results_file_train <- paste0("Contrast_", prop ,"_train.csv")
  results_file_test <- paste0("Contrast_", prop ,"_test.csv")
  results_file_inf <- paste0("Contrast_", prop ,"_inf.csv")

  fwrite(train_contrast_res,
        file.path(setup$output_directory, results_file_train))
  fwrite(test_contrast_res,
        file.path(setup$output_directory, results_file_test))
  fwrite(inf_contrast_res,
        file.path(setup$output_directory, results_file_inf))

  # Extract logFC for correlation analysis
  train_log_FC <- train_contrast_res |>
    select(coef, Feature, logFC) |>
    rename(logFC_train = logFC)

  test_log_FC <- test_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_test = logFC)

  inf_log_FC <- inf_contrast_res |>
    select(coef, Feature, logFC) |>
    rename(logFC_inf = logFC)

  # Merge all three data frames (aligned by coef and feature)
  merged_log_FC <- train_log_FC |>
    inner_join(test_log_FC, by = c("coef", "Feature")) |>
    inner_join(inf_log_FC, by = c("coef", "Feature"))

  # Correlation
  pearson <- head(compute_correlation(data = merged_log_FC, method = "pearson"), 2)
  spearman <- head(compute_correlation(data = merged_log_FC, method = "spearman"), 2)

  # Make metric dataframe
  metric_df <- inner_join(pearson, spearman, by = c("V1", "V2")) |>
    mutate(
      Seed = setup$seed,
      Status = ifelse(V1 == "logFC_test", "Test", "Inference"),
      Ancestry = toupper(setup$classification$infer_ancestry),
      Prediction = ifelse(V1 == "logFC_test", "Subset", "Ancestry"),
      n_train_ancestry = train_n,
      n_inf_ancestry = inf_n,
      proportion = gsub("\\_", ".", prop)
    )
  
  # Add 'metric_df' to list to combine later
  robustness_results <- append(robustness_results, list(metric_df))
}

# Combine into one big dataframe
final_metric_df <- bind_rows(robustness_results)

# Save metric Dataframe
fwrite(final_metric_df, file.path(setup$output_directory, "Metric_dge.csv"))

# Print statement to switch back to python script
print("Switching back to Python.")

