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
    library(fgsea)
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
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

} else {
  # Dev settings if no command-line argument provided
  yaml_file <- "job_settings.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}

# Save the results of the analysis: 'Interactions'
# 'Vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out = "data/combined_runs"
# Construct 'path_to_save_location'
tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
analysis_name = "Interactions"
path_to_save_location <- file.path(vscratch_dir_out, comparison, analysis_name)

# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Load the data
adata <- read_h5ad(setup$data_path)

# Define classification task
data <- adata[adata$obs[[setup$classification$output_column]]
              %in% setup$classification$comparison]

# The reference in the 'Interaction' analysis should be the same as in 'cross_ancestry'
# Setting intercept of ancestry: EUR
# Setting intercept of comparison: Based on occurance in 'setup$classification$comparison'
# 1. Transform settings into R useable form
output_column <- setup$classification$output_column
ancestry_column <- setup$classification$ancestry_column
train_ancestry <- setup$classification$train_ancestry

# Relevel 'ancestry_column' to be 'train_ancestry' as the intercept
# 1. Get possible values for ancestries
# 2. Relevel them so 'train_ancestry' is first level
# 3. Update in the actual data frame
possible_ancestry <- unique(data$obs[[ancestry_column]])
releveled_ancestry <- factor(possible_ancestry, 
                             levels = c(train_ancestry, setdiff(possible_ancestry, train_ancestry))
                            )
levels_ancestry <- levels(releveled_ancestry)

# Update 'ancestry_column' in the actual data
data$obs <- data$obs |>
  mutate(!!ancestry_column := factor(.data[[ancestry_column]], 
                                     levels = levels_ancestry))

# Relevel 'output_column' based on the occurance in the settings
# This is consitent with 'cross_ancestry' 
# 1. Get comparison from settings
# 2. Level based of occurance in settings
# 3. Update in the actual dataframe
possible_comparison <- unique(data$obs[[output_column]])
releveld_comparison <- factor(setup$classification$comparison)
levels_comparison <- levels(releveld_comparison)

# Update 'output_column' in the actual data
data$obs <- data$obs |>
  mutate(!!output_column := factor(.data[[output_column]], 
                                   levels = levels_comparison))

# Store baseline in case 
# 1. Find all columns in the data that are a factor
# 2. Construct the baseline 
reference_levels <- sapply(data$obs, function(col) {
  if (is.factor(col)) levels(col)[1] else NA
})

baseline <- paste(reference_levels[[ancestry_column]],
                  reference_levels[[output_column]],
                  sep = ":")

# Assertertion: 
# After a lot of manipulation check order of data
stopifnot(all(colnames(t(data$X)) ==
                row.names(data$obs[colnames(t(data$X)), ])))

# Create design matrix
# 1. Dynamically construct the formula
# 2. Create the design matrix
# 3. Rename colnames
formula <- as.formula(paste0("~", ancestry_column, " * ", output_column))
design <- model.matrix(formula, data = data$obs)

colnames(design) <- colnames(design) |>
  gsub(pattern = ancestry_column, replacement = "", colnames(design)) |> 
  gsub(pattern = output_column, replacement = "", colnames(design))

# TODO - Create DGEList object
# TODO - CalcNormFactors
# DGEList 
# 1. Calculate normfactors to make samples comparable
# 2. Filter lowly expressed genes
# dge <-  DGEList(counts = t(data$X))
# dge <- calcNormFactors(dge)

# Filter genes by expression
# TODO - Stronger filter maybe
keeper_genes <- filterByExpr(t(data$X), 
                             design = design,
                             min.count = 100,        # Minimum total counts across all samples
                             # min.total.count = 200, # Minimum total counts across all samples
                             # min.prop = 0.1,        # Require at least 10% of samples with non-zero counts
                             # min.cpm = 1,           # Minimum CPM of 1
                             # min.n = 3
                             )
filtered_data <- data[, keeper_genes]

# Limma workflow
# 1. Normalization 
data_norm <- voom(t(filtered_data$X), design, plot = FALSE)
# 2. Fit the model (including hypothesis testing)
# 3. Ebayes
# 4. Extract results
limma_fit <- lmFit(data_norm, design = design)
limma_fit <- eBayes(limma_fit)
limma_fit_res <- extract_results(limma_fit)
# Save the results
fwrite(limma_fit_res, file.path(path_to_save_location, "Interactions.csv"))

# Volcano plot 
# Showcase interesting interactions
volcano_plot <- limma_fit_res |>
  filter(str_detect(coef, ":")) |>
  mutate(coef = toupper(str_remove(coef, ":.*"))) |>
  ggplot(
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > 2)
    )
  ) +
  geom_point(alpha = 0.8, size = 3) +
  facet_grid(cols = vars(coef)) + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "grey"),
    name = 'Above thresholds'
  ) +
  labs(
    title = "Differences in genotype-phenotype relationship"
  )

# Save 'volcano_plot'
ggsave(file.path(path_to_save_location, "Volcano_plot.pdf"), 
       plot = volcano_plot, height = 5, width = 8)


# Find ancestry-specific enricht terms 
# Databases for gsea
db <- read_enrichR_database("data/geneset-libraries/MSigDB_Hallmark_2020.txt")
# 1. Filter for interaction effects
# 2. Split 'interactions' by ancestry
# 3. Run FGSEA
interactions <- limma_fit_res |> filter(str_detect(coef, ":"))
ancestry_specific_interactions <- interactions |> group_split(coef) 

# Function to run FGSEA
# Rank based on logFC
run_fgsea <- function(split_df, gene_sets) {
  # Create ranked list of genes
  ranked_genes <- split_df |>
    arrange(desc(logFC)) |>  # Order by logFC
    distinct(Feature, .keep_all = TRUE) |>  # Ensure unique genes
    pull(logFC, Feature)  # Pull named vector: logFC values with gene names
  
  # Perform FGSEA
  fgsea_results <- fgsea(pathways = gene_sets,
                         stats = ranked_genes
                         ) 
  
  # Return results as a tibble
  return(as_tibble(fgsea_results))
}

# Apply FGSEA to each split and store results
fgsea_results_list <- lapply(ancestry_specific_interactions, function(split) {
  coef_name <- unique(split$coef) # Name of the group
  results <- run_fgsea(split, db)
  results <- results |> mutate(coef = coef_name) # Add group info
  return(results)
})

# Bind list object into dataframe
fgsea_results <- bind_rows(fgsea_results_list)

# Visualize 
fgsea_results |>
  mutate(coef = toupper(str_remove(coef, ":.*"))) |>
  ggplot(
    aes(
      x = NES,
      y = pathway,
      size = -log10(padj)
    )
  ) + 
  geom_point() +
  facet_grid(
    cols = vars(coef)
  )











# Ranked gene list based on t-value
ranked_interactions <- setNames(ancestry_specific_interactiona[[1]]$t, ancestry_specific_interactiona[[1]]$Feature)
# Results
res <- fgsea(
  pathways = db,
  stats = ranked_interactions)



# Open pdf file
pdf(file.path(output_path, "Voom_plot.pdf"))
par(mfrow=c(1,2))
# Normalization (logCPM)
data_norm <- voom(t(data$X), design, plot = TRUE)
# Fit the model
limma_fit <- lmFit(data_norm, design = design)
# Ebayes
limma_fit <- eBayes(limma_fit)
# Plot
plotSA(limma_fit, main="Final model: Mean-variance trend")
# Close the PDF device
dev.off()

# Extract results
limma_res <- list()
for (x in colnames(coef(limma_fit))){
  # Extract for each coefficient all genes
  limma_res[[x]] <- topTable(limma_fit, coef = x, number = Inf) |>
     rownames_to_column("Feature")
}
limma_res <- bind_rows(limma_res, .id = "coef")
limma_res <- filter(limma_res, coef != "(Intercept)")

# Save the results from statistics
fwrite(limma_res, 
       file.path(output_path, "Interaction_results.csv"))

# Filter for interaction terms
limma_res_interactions <- filter(limma_res, str_detect(coef, ":"))

# limma_res_sig |>
#   filter(str_detect(coef, ":")) |>
#   group_by(coef) |>
#   count()

# limma_res_sig |>
#   filter(str_detect(coef, ":")) |>
#   filter(logFC < 0) |>
#   mutate(coef = str_replace(coef, ":.*$", "")) |>
#   filter(coef == 'afr') |>
#   pull(Feature)



# Plot them as a volcano plot
volcano_plot <- limma_res_interactions |>
  ggplot(
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > 2)
    )
  ) +
  geom_point(alpha = 0.8, size = 3) +
  facet_grid(cols = vars(coef)) + 
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "grey"),
    name = 'Above thresholds'
  ) 

# Save the volcano plot
ggsave(file.path(output_path, "Volcano_plot.pdf"),
       plot = volcano_plot, height = 6, width = 12)

# Pvalue distribution
p_value_dist <- limma_res |>
  ggplot(aes(x = P.Value,
  fill=factor(floor(AveExpr)))) +
  geom_histogram() +
  facet_grid(rows = vars(coef))

# Save the distribution plot
ggsave(file.path(output_path, 
       "P_value_plot.pdf"),
       plot = p_value_dist, height = 10, width = 6)

# Plot for interactions
p_value_dist <- limma_res_interactions |>
  ggplot(aes(x = P.Value,
  fill=factor(floor(AveExpr)))) +
  geom_histogram() +
  facet_grid(rows = vars(coef))

# Save the distribution plot
ggsave(file.path(output_path,
       "P_value_interactions_plot.pdf"),
       plot = p_value_dist, height = 10, width = 6)

# Filter significant interactions
limma_res_sig <- limma_res |>
  as_tibble() |>
  filter(adj.P.Val < 0.05)

# Plot feature with biggest interaction effect
feature <- limma_res_sig |>
  filter(str_detect(coef, ":")) |>
  group_by(coef) |>
  slice_max(logFC) |>
  pull(Feature)

# Get statistic
# annotation <- limma_res_sig |>
#   filter(str_detect(coef, ":")) |>
#   group_by(coef) |>
#   slice_max(logFC) |>
#   select(coef, Feature, logFC , adj.P.Val) |>
#   mutate(coef = gsub(":Breast_Invasive_Lobular_Carcinoma", "", coef)) |>
#   mutate({{ancestry_column}} := coef)

# Plot this one feature
feature_expression <- as_tibble(t(data_norm$E[feature,]), rownames = 'sampleId')
meta_data <- as_tibble(data$obs, rownames = 'sampleId')
feature_meta <- merge(meta_data, feature_expression, on = 'sampleId')
# Pivot longer
feature_meta <- feature_meta |>
  pivot_longer(cols = all_of(feature),
               names_to = "Feature",
               values_to = "Value")

# Plot
single_feature_plot <- feature_meta |>
  ggplot(
    aes(
      x = get(output_column),
      y = Value
    )
  ) +
  geom_boxplot() +
  facet_grid(cols = vars(get(ancestry_column)),
             rows = vars(Feature)) +
  xlab("Comparsion") + 
  ylab("Normalized Expression") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(output_path,
       "Single_features.pdf"),
       plot = p_value_dist, height = 10, width = 6)


# Visualize multiple features
# Get up reg. features
goi_up <- limma_res_sig |>
  filter(str_detect(coef, ":"), 
         adj.P.Val < 0.05) |>
  group_by(coef) |>
  filter(logFC > 0) |>
  slice_max(logFC, n = 5) |>
  pull(Feature)

# Get down reg. features
goi_down <- limma_res_sig |>
  filter(str_detect(coef, ":"), 
         adj.P.Val < 0.05) |>
  group_by(coef) |>
  filter(logFC < 0) |>
  slice_min(logFC, n = 5) |>
  pull(Feature)

# Statistic results
# Plot up reg. features
stat_up <- limma_res_interactions |>
  filter(Feature %in% goi_up) |>
  ggplot(
    aes(
      y = Feature,
      x = str_replace(coef, ":.*$", ""),
      color = logFC,
      size = -log10(adj.P.Val))
  ) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  xlab("Comparison") +
  ylab("Feature") +
  ggtitle("Top upreg. interactions") + 
  theme(axis.text.x = element_text(angle = 90))

# Expression
exp_list <- list()
for(gg in goi_up){
  exp_list[[gg]] <- data$obs |>
    mutate(E=scale(data_norm$E[gg,])) |>
    rownames_to_column("sampleId") |>
    remove_rownames()
}
# Bind to dataframe
expression_up <- bind_rows(exp_list, .id = "Feature")

# Plot expression heatmap
heatmap_up <- expression_up |>
  group_by(cancer_type_detailed, genetic_ancestry, Feature) |>
  summarise(across(where(is.numeric), mean)) |>
  ggplot(
    aes(
      x = genetic_ancestry,
      y = Feature,
      fill = E
    )
  ) +
  geom_tile() +
  facet_grid(cols = vars(cancer_type_detailed)) +
  scale_fill_gradient2(low="blue", high="red") +
  ggtitle("Average expression per ancestry") 

# Patchwork (combine stats results and expression)
upregulated_plot <- heatmap_up + stat_up
# Save the patchwork plot
ggsave(file.path(output_path, 
       "Upregulated_interactions.pdf"),
       plot = upregulated_plot, height = 5, width = 10)

# Plot down reg. features
stat_down <- limma_res_interactions |>
  filter(Feature %in% goi_down) |>
  ggplot(
    aes(
      y = Feature,
      x = str_replace(coef, ":.*$", ""),
      color = logFC,
      size = -log10(adj.P.Val))
  ) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  xlab("Comparison") +
  ylab("Feature") +
  ggtitle("Top downreg. interactions") + 
  theme(axis.text.x = element_text(angle = 90))

# Expression
exp_list <- list()
for(gg in goi_down){
  exp_list[[gg]] <- data$obs |>
    mutate(E=scale(data_norm$E[gg,])) |>
    rownames_to_column("sampleId") |>
    remove_rownames()
}
# Bind to dataframe
expression_down <- bind_rows(exp_list, .id = "Feature")

# Plot expression heatmap
heatmap_down <- expression_down |>
  group_by(cancer_type_detailed, genetic_ancestry, Feature) |>
  summarise(across(where(is.numeric), mean)) |>
  ggplot(
    aes(
      x = genetic_ancestry,
      y = Feature,
      fill = E
    )
  ) +
  geom_tile() +
  facet_grid(cols = vars(cancer_type_detailed)) +
  scale_fill_gradient2(low="blue", high="red") +
  ggtitle("Average expression per ancestry") 

# Patchwork (combine stats results and expression)
downregulated_plot <- heatmap_down + stat_down
# Save the patchwork plot
ggsave(file.path(output_path, 
       "Downregulated_interactions.pdf"),
       plot = downregulated_plot, height = 5, width = 10)









