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
robustness_results <- list()
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
  stopifnot(all(colnames(t(train_data$X)) ==
                  row.names(train_data$obs[colnames(t(train_data$X)), ])))

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

  # Create design matrix
  train_design <- create_design(setup$classification$output_column,
                                meta = train_data$obs,
                                covariate = covariate_list)                            
  test_design <- create_design(setup$classification$output_column,
                              meta = test_data$obs,
                              covariate = covariate_list)
  inf_design <- create_design(setup$classification$output_column,
                              meta = inf_data$obs,
                              covariate = covariate_list)

  # Filter genes by expression (Use the train set)
  # This will lead to different features across proportions
  # TODO - This will lead to different features across proportions
  keeper_genes <- filterByExpr(t(train_data$X), train_design)
  # Subset features
  train_filtered <- train_data[, keeper_genes]
  test_filtered <- test_data[, keeper_genes]
  inf_filtered <- inf_data[, keeper_genes]

  # Save features (for use in ML)
  # Features are saved per 'prop'
  feature_file <- paste0("Features_", prop_, ".yml")
  write_yaml(train_filtered$var_names, file.path(setup$output_directory, feature_file))

  # Assertion: Check array
  # Same genes in all sets
  stopifnot(all(train_data$var_names == test_data$var_names))
  # Dimensions
  stopifnot(table(keeper_genes)[2] == dim(train_filtered$X)[2])

  # Limma workflow
  print("Start differential gene expression analysis.")
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

  # Assertion: Check if all genes have been tested
  stopifnot(all(dim(train_mean_res)[1] / length(unique(train_mean_res$coef)) ==
                table(keeper_genes)[2]))

  # Save results
  fwrite(train_mean_res, file.path(setup$output_directory, paste0("Means_", prop_, "_train.csv")))
  fwrite(test_mean_res, file.path(setup$output_directory, paste0("Means_", prop_, "_test.csv")))
  fwrite(inf_mean_res, file.path(setup$output_directory, paste0("Means_", prop_, "_inf.csv")))

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
  results_file_train <- paste0("Contrast_", prop_ ,"_train.csv")
  results_file_test <- paste0("Contrast_", prop_ ,"_test.csv")
  results_file_inf <- paste0("Contrast_", prop_ ,"_inf.csv")

  fwrite(train_contrast_res,
        file.path(setup$output_directory, results_file_train))
  fwrite(test_contrast_res,
        file.path(setup$output_directory, results_file_test))
  fwrite(inf_contrast_res,
        file.path(setup$output_directory, results_file_inf))

  # Filter coef that are in comparison (only interested in those)
  train_contrast_res <- train_contrast_res |> filter(coef %in% setup$classification$comparison)
  test_contrast_res <- test_contrast_res |> filter(coef %in% setup$classification$comparison)
  inf_contrast_res <- inf_contrast_res |> filter(coef %in% setup$classification$comparison)

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
  
  # Add 'metric_df' to list to combine later
  robustness_results <- append(robustness_results, list(metric_df))
}

# Combine into one big dataframe
final_metric_df <- bind_rows(robustness_results)

# Save metric Dataframe
fwrite(final_metric_df, file.path(setup$output_directory, "Metric_dge.csv"))

# Print statement to switch back to python script
print("Switching back to Python.")

