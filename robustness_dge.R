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
# Set seed (because why not)
set.seed(setup$seed)

# Load data
adata <- read_h5ad(setup$data_path)

# Iterate over all proportions doing DGE analysis
robustness_metric <- list()
robustness_logFC <- list()
for (prop in setup$proportion){

  print(paste("Dealing with", prop, "% of data."))
  # Replace . with _ in 'prop'
  prop_ <- gsub("\\.", "_", prop)

  # Create file names
  observation_file_train <- paste0("Obs_proportion_", prop_, "_train.yml")
  observation_file_test <- paste0("Obs_proportion_", prop_, "_test.yml")
  observation_file_inf <- paste0("Obs_proportion_", prop_, "_inf.yml")

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
  stopifnot(all(colnames(t(train_data$X)) == row.names(train_data$obs[colnames(t(train_data$X)), ])))

  # Covariate:
  # 1. Extract the covariate, if it exists and has a meaningful value
  covariate_list <- get_covariate(setup)
  # 2. Check if there are replicates and at least two levels for the covariate 
  if (!is.null(covariate_list)){
    # Check covariate type
    print("Checking covariates.")
    covariate_types <- classify_covariates(train_data$obs, covariate_list)
    # Check conditions "Discrete"
    check_covariate_conditions(train_data$obs, covariates = covariate_types$discrete)
    check_covariate_conditions(test_data$obs, covariates = covariate_types$discrete)
    check_covariate_conditions(inf_data$obs, covariates = covariate_types$discrete)

    # Visualize
    if (setup$visual_val){
    # Visualize all discrete columns
    to_visualize_columns  <- c(setup$classification$output_column, covariate_types$discrete)
    # Datasets to visualize
    datasets <- list(
    Train = train_data$obs,
    Test = test_data$obs,
    `Inf` = inf_data$obs
      )
    # Function that creates the plots
    patchwork_plot <- create_stratification_plot(
      datasets = datasets,
      to_visualize_columns = to_visualize_columns
      )
    # Save
    ggsave(
      filename = paste0("Validation_", prop_, "_stratification.pdf"),
      plot = patchwork_plot,
      path = setup$output_directory,
      width = 5, height = 5
      )
    }
  } 

  # Create design (GLM)
  output_column <- setup$classification$output_column
  train_design <- create_design(output_column, train_data$obs, covariate = covariate_list)                            
  test_design <- create_design(output_column, test_data$obs, covariate = covariate_list)
  inf_design <- create_design(output_column, inf_data$obs, covariate = covariate_list)

  # Filter genes
  if (setup$filter_features & setup$data_type == "expression"){
    print("Filter features")
    # Filter by expression
    keeper_genes <- filterByExpr(t(train_data$X), design = train_design)
    # Subset features
    train_filtered <- train_data[, keeper_genes]
    test_filtered <- test_data[, keeper_genes]
    inf_filtered <- inf_data[, keeper_genes]
  } else{
    train_filtered <- train_data
    test_filtered <- test_data
    inf_filtered <- inf_data
  }
  # Save features (for use in ML)
  output_directory <- setup$output_directory
  file_name <- paste0("Features_", prop_, ".yml")
  write_yaml(train_filtered$var_names, file.path(output_directory, file_name))

  # Limma workflow
  print("Start differential gene expression analysis.")
  # Select normalization method
  data_type <- setup$data_type
  dge_normalization <- setup$dge_normalization
  normalization_method <- normalization_methods[[data_type]][[dge_normalization]]
  
  # Transpose (rows = Genes, cols = Samples)
  train_filtered = t(train_filtered$X)
  test_filtered = t(test_filtered$X)
  inf_filtered = t(inf_filtered$X)

  # Normalization
  train_norm <- normalization_method(train_filtered, train_design)
  test_norm <- normalization_method(test_filtered, test_design)
  inf_norm <- normalization_method(inf_filtered, inf_design)

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

  # Save results
  output_directory <- setup$output_directory
  file_name <- paste0("Limma_means_", prop_, "_train.csv")
  fwrite(train_mean_res, file.path(output_directory, file_name))

  file_name <- paste0("Limma_means_", prop_, "_test.csv")
  fwrite(test_mean_res, file.path(output_directory, file_name))

  file_name <- paste0("Limma_means_", prop_, "_inf.csv")
  fwrite(inf_mean_res, file.path(output_directory, file_name))

  # Create contrast matrix 
  print("Fit contrast matrix.")
  # Used for hypothesis testing between groups
  comparison <- setup$classification$comparison
  contrast_matrix_train <- create_contrast(colnames(train_design), conditions = comparison)
  contrast_matrix_test <- create_contrast(colnames(test_design), conditions = comparison)
  contrast_matrix_inf <- create_contrast(colnames(inf_design), conditions = comparison)

  # Fit contrast
  train_limma_fit_contrast <- contrasts.fit(train_limma_fit, contrast_matrix_train)
  test_limma_fit_contrast <- contrasts.fit(test_limma_fit, contrast_matrix_test)
  inf_limma_fit_contrast <- contrasts.fit(inf_limma_fit, contrast_matrix_inf)

  # Ebayes
  train_limma_fit_contrast <- eBayes(train_limma_fit_contrast)
  test_limma_fit_contrast <- eBayes(test_limma_fit_contrast)
  inf_limma_fit_contrast <- eBayes(inf_limma_fit_contrast)

  # Results contrast
  train_contrast_res <- extract_results(train_limma_fit_contrast)
  test_contrast_res <- extract_results(test_limma_fit_contrast)
  inf_contrast_res <- extract_results(inf_limma_fit_contrast)

  # Save contrast results
  output_directory <- setup$output_directory
  file_name <- paste0("Limma_contrast_", prop_ ,"_train.csv")
  fwrite(train_contrast_res, file.path(output_directory, file_name))

  file_name <- paste0("Limma_contrast_", prop_ ,"_test.csv")
  fwrite(test_contrast_res, file.path(output_directory, file_name))

  file_name <- paste0("Limma_contrast_", prop_ ,"_inf.csv")
  fwrite(inf_contrast_res, file.path(output_directory, file_name))

  # Filter coef that are in comparison (only interested in those)
  train_contrast_res <- filter(train_contrast_res, coef %in% setup$classification$comparison)
  test_contrast_res <- filter(test_contrast_res, coef %in% setup$classification$comparison)
  inf_contrast_res <- filter(inf_contrast_res,coef %in% setup$classification$comparison)

  # Foramatting results
  # Raw logFC
  train_log_FC <- train_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Train",
    Ancestry = toupper(setup$classification$train_ancestry),
  )

  test_log_FC <- test_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Test",
    Ancestry = toupper(setup$classification$train_ancestry),
  )

  inf_log_FC <- inf_contrast_res |>
  select(coef, Feature, logFC) |>
  mutate(
    Status = "Inference",
    Ancestry = toupper(setup$classification$infer_ancestry)
  )

  # Combine across sets
  raw_logFC <- bind_rows(train_log_FC, test_log_FC, inf_log_FC) |>
  mutate(
    Seed = setup$seed,
    n_train_ancestry = train_n,
    n_inf_ancestry = inf_n,
    Proportion = prop
  )

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
      Proportion = prop
    )
  
  # Combine across proportions
  robustness_metric <- append(robustness_metric, list(metric_df))
  robustness_logFC <- append(robustness_logFC, list(raw_logFC))
}

# Combine into one big dataframe
robustness_metric <- bind_rows(robustness_metric)
robustness_logFC <- bind_rows(robustness_logFC)

# Save 
output_directory <- setup$output_directory
fwrite(robustness_metric, file.path(output_directory, "Contrast_metric_dge.csv"))
fwrite(robustness_logFC, file.path(output_directory, "LogFCs.csv"))

# Print statement to switch back to python script
print("Switching back to Python.")

