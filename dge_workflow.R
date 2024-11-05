# Remove start up messages
suppressPackageStartupMessages(
  {
  # Standard libraries
  library(tidyverse)
  library(data.table)
  library(yaml)
  library(anndata)
  # Statistics library
  library(edgeR)
  # Visualization
  library(httpgd)
  }
)

# Custom functions
source('r_utils.R')

# Give permission to execute
# TODO - Q.Wolfgang: Is there a way to give the script execution permission?
system('chmod 777 dge_workflow.R')

# Here starts the script
print('Envoking R.')
print('Starting differential gene expression analysis.')

# Load the data (In not anndata format)
# meta <- fread('data/downloads/tcga_csv/tcga_studies_meta_intersection_adjusted_withnewcolumns.csv')[,-1]
# data <- fread('data/downloads/tcga_csv/seqV2_merged_RSEM.csv')[,-1]

# Meta data needs to be modified to replace all ' ' with '_'
# charvars <- sapply(meta, is.character)
# meta[,
#      (names(meta)[charvars]) := lapply(.SD, gsub, pat="[: ]", rep="_"),
#      .SDcols=charvars
# ]

# Load the settings (They are a command line input)
args = commandArgs(trailingOnly = TRUE)

# TODO - Check that it is a yaml file
setup <- yaml.load_file(args[1])

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

# Filter genes by expression (Use the train sets)
keeper_genes <- filterByExpr(t(train_data$X), train_design)
# Subset 
train_filtered <- train_data[, keeper_genes]
test_filtered <- test_data[, keeper_genes]
inf_filtered <- inf_data[, keeper_genes]

# Save feature names (For use in ML)
write_yaml(train_filtered$var_names, file.path(setup$output_directory, 'Features.yml'))

# Assertion: Check array
## Q.Wolfgang: How to check 3 list are equal ----
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
# TODO - check means model
# TODO - check contrast





print('Switching back to Python.')

# Create output directory
# make.dir(setup$output_directory, overwrite = setup$overwrite)

# Define training and inferred ancestry
# eur_meta <- meta[consensus_ancestry == setup$classification$train_ancestry, ]
# inf_meta <- meta[consensus_ancestry == setup$classification$infer_ancestry, ]
# eur_data <- adata[adata$obs['genetic_ancestry'] == setup$classification$train_ancestry]
# inf_data <- adata[adata$obs['genetic_ancestry'] == setup$classification$inf_ancestry]

# Define classification 
## eval(as.name(...)) was used to have input from setup file
# eur_meta <- eur_meta[eval(as.name(setup$classification$output_column)) 
                     # %in% setup$classification$comparison, ]
# inf_meta <- inf_meta[eval(as.name(setup$classification$output_column)) 
                     # %in% setup$classification$comparison, ]


# Assertion: Check for enough samples per class
# stopifnot(all(eur_meta[, .N, by = cancer_type_detailed] > 25))

# Set the seed
# set.seed(setup$seed)


# Normalization: in R use edgeR (filter by low expression) and voom
# Switch to tibble (Convenience)
# eur_meta <- as_tibble(eur_meta) |> 
# column_to_rownames(var = 'sampleId')
# inf_meta <- as_tibble(inf_meta) |> 
#  column_to_rownames('sampleId')
# Transpose sample id = columns, genes = index (for DGE workflow)
# data <- as_tibble(data) |> 
#   column_to_rownames(var = 'sampleId') |> 
#   t() 
# Filter data based on sample id 
# counts <- data[, rownames(eur_meta)]

# Assertion: Check arrays
# Dimensions
# stopifnot(dim(counts)[2] == dim(eur_meta)[1])
# Order
# stopifnot(all(colnames(counts) == row.names(eur_meta[colnames(counts), ])))

# Create the design matrix and contrast matrix
# Design
# design <- model.matrix(~0 + eval(as.name(setup$classification$output_column)), eur_meta) 
# Rename the columns 
# repl_str <- 'eval(as.name(setup$classification$output_column))'
# colnames(design) <- gsub(repl_str, '', colnames(design), fixed=TRUE)

# Contrast
# Compare the average of one condition to the average of all other
# OVR
# contrast.matrix <- create_contrast(colnames(design))

# Filter genes by expression
# filtered_genes <- filterByExpr(counts, design)
# filtered_counts <- counts[filtered_genes, ]

# Assertion: Check arrays
# Dimensions
# stopifnot(table(filtered_genes)[2] == dim(filtered_counts)[1])
# Order
# stopifnot(all(colnames(filtered_counts) == row.names(eur_meta[colnames(filtered_counts), ])))

# Normalization
# norm_counts <- voom(filtered_counts, design, plot = TRUE)





