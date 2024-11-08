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
  use_condaenv('/opt/conda/envs/ancestry/bin/python')
  # Statistics library
  library(edgeR)
  # Visualization
  library(patchwork)
  library(httpgd)
  }
)

# Custom functions
source('r_utils.R')

# Here starts the script
print('Envoking R.')
# Load the settings (They are a command line input)
args = commandArgs(trailingOnly = TRUE)
# Check that it is a yaml file and exists
is_yml_file(args[1])
setup <- yaml.load_file(args[1])


# Data loading for debugging
setup <- yaml.load_file('/home/people/dkatzlberger/Project_Y/dev_settings.yml')

print('Start differential gene expression analysis.')

# Load data
adata <- read_h5ad(setup$data_path)

train_idx <- yaml.load_file(file.path(setup$output_directory, 'Obs_train.yml'))
test_idx <- yaml.load_file(file.path(setup$output_directory, 'Obs_test.yml'))
inf_idx <- yaml.load_file(file.path(setup$output_directory, 'Obs_inf.yml'))

# Subset the data by the indexes create in python
train_data <- adata[adata$obs_names %in% train_idx, ]
test_data <- adata[adata$obs_names %in% test_idx, ]
inf_data <- adata[adata$obs_names %in% inf_idx, ]

# Assertion: Check arrays
# Dimensions
stopifnot(length(train_idx) == length(train_data$obs_names))
# Order
stopifnot(all(colnames(t(train_data$X)) == row.names(train_data$obs[colnames(t(train_data$X)), ])))

# Create design matrix
train_design <- create_design(setup$classification$output_column, meta = train_data$obs)
test_design <- create_design(setup$classification$output_column, meta = test_data$obs)
inf_design <- create_design(setup$classification$output_column, meta = inf_data$obs)


# Create contrast matrix (Only needs to be created once)
contrast_matrix <- create_contrast(colnames(train_design))
print(contrast_matrix)

# Filter genes by expression (Use the train sets)
keeper_genes <- filterByExpr(t(train_data$X), train_design)
# Subset 
train_filtered <- train_data[, keeper_genes]
test_filtered <- test_data[, keeper_genes]
inf_filtered <- inf_data[, keeper_genes]

# Save feature names (For use in ML)
write_yaml(train_filtered$var_names, file.path(setup$output_directory, 'Features.yml'))

# Assertion: Check array
# Same genes in all sets
stopifnot(all(train_data$var_names == test_data$var_names))
# Dimensions
stopifnot(table(keeper_genes)[2] == dim(train_filtered$X)[2])

# Normalization (Voom)
train_norm <- voom(t(train_filtered$X), train_design, plot = FALSE)
test_norm <- voom(t(test_filtered$X), test_design, plot = FALSE)
inf_norm <- voom(t(inf_filtered$X), inf_design, plot = FALSE)

# Fit the model (Means model)
train_limmaFit <- lmFit(train_norm, design=train_design)
test_limmaFit <- lmFit(test_norm, design=test_design)
inf_limmaFit <- lmFit(inf_norm, design=inf_design)

# Ebayes
train_limmaFit <- eBayes(train_limmaFit)
test_limmaFit <- eBayes(test_limmaFit)
inf_limmaFit <- eBayes(inf_limmaFit)

# Results means model
train_mean_res <- extract_results(train_limmaFit)
test_mean_res <- extract_results(test_limmaFit)
inf_mean_res <- extract_results(inf_limmaFit)

# Assertion: Check if all genes have been tested
stopifnot(all(dim(train_mean_res)[1]/length(unique(train_mean_res$coef)) == table(keeper_genes)[2]))

# Save results
fwrite(train_mean_res, file.path(setup$output_directory, 'Means_train.csv'))
fwrite(test_mean_res, file.path(setup$output_directory, 'Means_test.csv'))
fwrite(inf_mean_res, file.path(setup$output_directory, 'Means_inf.csv'))

# Fit contrast
train_limmaFit_contrast <- contrasts.fit(train_limmaFit, contrast_matrix)
test_limmaFit_contrast <- contrasts.fit(test_limmaFit, contrast_matrix)
inf_limmaFit_contrast <- contrasts.fit(inf_limmaFit, contrast_matrix)

# Ebayes
train_limmaFit_contrast <- eBayes(train_limmaFit_contrast)
test_limmaFit_contrast <- eBayes(test_limmaFit_contrast)
inf_limmaFit_contrast <- eBayes(inf_limmaFit_contrast)

# Results contrast
train_contrast_res <- extract_results(train_limmaFit_contrast)
test_contrast_res <- extract_results(test_limmaFit_contrast)
inf_contrast_res <- extract_results(inf_limmaFit_contrast)

# Save results
fwrite(train_contrast_res, file.path(setup$output_directory, 'Contrast_train.csv'))
fwrite(test_contrast_res, file.path(setup$output_directory, 'Contrast_test.csv'))
fwrite(inf_contrast_res, file.path(setup$output_directory, 'Contrast_inf.csv'))

# Visualization
# Check if logFC of models align with signal in the data
if(setup$visual_val){
# Get genes for validation
top <- train_mean_res |> 
  slice_max(logFC, n = 5) |> 
  pull(Gene) |> 
  unique()
  
low <- train_mean_res |> 
  slice_min(logFC, n = 5) |> 
  pull(Gene) |> 
  unique()

goi <- c(top, low)

# Validation of the training 
t <- visual_validation(
  meta = train_data$obs,
  signal = train_norm$E,
  mean_stats = train_mean_res,
  contrast_stats = train_contrast_res,
  goi = goi,
  data_output_column = setup$classification$output_column
)

i <- visual_validation(
  meta = inf_data$obs,
  signal = inf_norm$E,
  mean_stats = inf_mean_res,
  contrast_stats = inf_contrast_res,
  goi = goi,
  data_output_column = setup$classification$output_column
)

# Save the pictures
  ggsave(file.path(setup$output_directory, 'Validation_stats_train.pdf'), plot = t, height = 6, width = 12)
  ggsave(file.path(setup$output_directory, 'Validation_stats_inf.pdf'), plot = i, height = 6, width = 12)
}

# Print statement to switch back to python script
print('Switching back to Python.')

# x <- function(s){
# train_contrast_res |> 
#   select(c(all_of(s), Gene)) |> 
#   group_by(c(s)) |> 
#   pivot_longer(
#     cols = !c(s),
#     names_to = 'name',
#     values_to = 'Gene'
#   )
# }


# x('coef')

