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
  yaml_file <- "job_settings.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
}

# Save directory
# Is often modified using a script which manipulates path to compute many seeds at once
path_to_save_location <- setup$output_directory
# Create directory if not exist inclusive parents
if (!dir.exists(directory)) {
  dir.create(directory, recursive = TRUE)
}

# This script relies on random sampling
# Hence a 'seed' for reproducability
seed <- setup$seed
set.seed(setup$seed)

# Load the data (Anndata format)
adata <- read_h5ad(setup$data_path)

# Define classification task
data <- adata[adata$obs[[setup$classification$output_column]]
              %in% setup$classification$comparison]
# Filter only EUR based on ancestry column (column of used ancestry)
eur_data <- data[data$obs[[setup$classification$ancestry_column]] 
                 == setup$classification$train_ancestry]

# Randomly split EUR into two subsets (without stratification)
# 1. Sample random index to split the data frame (50 50)
# 2. Subset data based on these indexes 
split_index <- sample(1:nrow(eur_data), size = nrow(eur_data) / 2)
# a) The 'constant_split' is keept as ground_truth (no further splitting)
# b) The 'sampling_split' is used to randomly sample by size
#   i) The only prerequesite is that for each class there are at least two replicates
#      E.g. Cancer_1 = 2 observations, Cancer_2 = 2 observations
#      This is needed by limma eBayes
constant_split <- eur_data[split_index, ]
sampling_split <- eur_data[-split_index, ]

# Calculate logFC for the 'constant_split' (ground_truth)
# 1. Creating design matrix for the 'constant split'
# 2. Create contrast matrix for hyptheses testing (only done once)
design <- create_design(setup$classification$output_column, meta = constant_split$obs)
contrast_matrix <- create_contrast(colnames(design))

# Filter lowley expressed genes bases on 'constant_split'
keeper_genes <- filterByExpr(t(constant_split$X), design)
# Remove lowley expressed genes
# 1. For the 'constant_split'
# 2. For the 'sampling_split'
constant_split <- constant_split[, keeper_genes]
sampling_split <- sampling_split[, keeper_genes]

# TODO - DGEList object for 'constant_split'
# TODO - Add calcNormfactors

# Normalization and logCPM transformation 'constant_split'
constant_norm <- voom(t(constant_split$X), design, plot = FALSE)
# Limma
# 1. Fit the model (means model)
# 2. Ebayes
# 3. Extract results
constant_limma_fit <- lmFit(constant_norm, design = design)
constant_limma_fit <- eBayes(constant_limma_fit)
constant_means_res <- extract_results(constant_limma_fit)

# Hypothesis testing
# 1. Fit contrast
# 2. Ebayes
# 3. Extreact results
constant_limma_fit_contrast <- contrasts.fit(constant_limma_fit, contrast = contrast_matrix)
constant_limma_fit_contrast <- eBayes(constant_limma_fit_contrast)
constant_contrast_res <- extract_results(constant_limma_fit_contrast)

# From the Toptable only need logFCs are interesting
# Rename the logFC column to indicate they are from the 'constant_split'
constant_logFC <- constant_contrast_res |>
  select(coef, Feature, logFC) |>
  rename(logFC_constant = logFC) 

# Randomly sample the 'sampling_split'
# 1. Define sample size
# 2. A subset by size is randomly sample from 'sampling_split'

# The sample size is calculated based on total n of the 'sampling_split'
# Remember: 'Sampling_split' has the same n as 'constant_split'
props <- c(1.0, 0.8, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02)
sizes <- floor(props * nrow(sampling_split))

# 'Sample_by_size' only ensures that at least 2 replicates per class are present in each subset
# No stratification
subset_list <- sample_by_size(sampling_split, 
                              sizes = sizes,
                              class_column = setup$classification$output_column,
                              seed = setup$seed)

# Visualize each sampled subset 
sample_meta <- data.frame()
for (i in seq_along(subset_list)) {
  smpl <- subset_list[[i]]
  sample_name <- names(subset_list)[i]
  # Get meta and add column Sample
  meta <- smpl$obs |> 
    mutate(Sample = sample_name)
  # Combine
  sample_meta <- bind_rows(sample_meta, meta) 
}
# Add meta of 'stable_subset'
sample_meta <- 
  constant_split$obs |>
  mutate(Sample = "constant_observations") |>
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

# Save the pictures
ggsave(file.path(setup$output_directory, "Subsets.pdf"),
        plot = sample_distribution, height = 6, width = 12)

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

# Make metric dataframe with all important information 
metric_df <- inner_join(pearson, spearman, by = c("V1", "V2")) |>
  mutate(
    Ancestry = "EUR"
    Seed = setup$seed,
    Proportions = props,
    n_proportion = sizes,
    n_stable = stable_n_obs
  )

# Save metric Dataframe
fwrite(metric_df, file.path(setup$output_directory, "Metric_subsetting.csv"))