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

# Start with 'constant_split'
# 1. Load observations
constant_index <- yaml.load_file(file.path(setup$output_directory, "Obs_constant.yml"))
constant_split <- adata[adata$obs_names %in% constant_index, ]

# Calculate logFC for the 'constant_split' (ground_truth)
# 1. Creating design matrix for the 'constant split'
# 2. Create contrast matrix for hyptheses testing (only done once)
design <- create_design(setup$classification$output_column, meta = constant_split$obs)
contrast_matrix <- create_contrast(colnames(design))

# Filter lowley expressed genes bases on 'constant_split'
keeper_genes <- filterByExpr(t(constant_split$X), design)
constant_split <- constant_split[, keeper_genes]

# Save feature names (For use in ML)
write_yaml(constant_split$var_names,
           file.path(setup$output_directory, "Features.yml"))

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


# Load subsets
vscratch_dir <- setup$output_directory
match_pattern <- paste0("Obs_proportion_")

# Extracting all folders in the 'vscratch_dir' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir'
# 2. With 'match_pattern' extract matching folders
all_vscratch_dir <- list.files(vscratch_dir, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir, value = TRUE)

# Create a list with all subsets
# 1. Filters 'adata' for observations
# 2. Filters 'adata' for features
subset_list <- list()
for (file in match_vscratch_dir){
    index <- yaml.load_file(file)
    # Subset by index and filtered features
    subset <- adata[adata$obs_names %in% index, keeper_genes]
    # Append to list
    subset_list <- append(subset_list, subset)
}

# Limma analysis for all other subsets
# 1. Initalize a instance to store results
# 2. Iterate over all subsets in 'subset_list'
logFC_results <- NULL
for (i in seq_along(subset_list)){
  subset <- subset_list[[i]]
  # Check if observations in subset are unique
  stopifnot(anyDuplicated(subset$obs_names) == 0)

  # Information about the subset
  n_obs <- nrow(subset)

  # TODO - Create DGEList object
  # TODO - CalcNormFactors

  # Limma workflow
  subset_design <- create_design(setup$classification$output_column, meta = subset$obs)
  subset_norm <- voom(t(subset$X), subset_design, plot = FALSE)
  # 3. Fit the model (means model)
  # 4. Ebayes
  # 5. Extract results
  subset_limma_fit <- lmFit(subset_norm, design = subset_design)
  subset_limma_fit <- eBayes(subset_limma_fit)
  susbet_means_res <- extract_results(subset_limma_fit)

  # Hypothesis testing
  # 1. Fit contrast
  # 2. Ebayes
  # 3. Extract results
  subset_limma_fit_contrast <- contrasts.fit(subset_limma_fit, contrast_matrix)
  subset_limma_fit_contrast <- eBayes(subset_limma_fit_contrast)
  subset_contrast_res <- extract_results(subset_limma_fit_contrast)

  # Extract logFC for correlation analysis
  # Uses the subset name to distinguish 
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

# Correlation of all numeric columns (logFC)
# 1. Pearson
# 2. Spearman
# Only keep correlation with 'constant_split'
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
    n = str_extract(V1, "\\d+$")
  )

# Save metric Dataframe
fwrite(metric_df, file.path(setup$output_directory, "Metric_eur_subsetting_dge.csv"))


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
    axis.text.x = element_text(angle = 90),
    plot.caption = element_text(hjust = 0)
  ) 

# Save 'subset_plot'
ggsave(file.path(setup$output_directory, "Subsets.pdf"),
       plot = subset_plot, height = 6, width = 12)
}

# Print statement to switch back to python script
print("Switching back to Python.")
