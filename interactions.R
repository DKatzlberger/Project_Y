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
base_dir = "data/combined_runs"
tag <- setup$tag
comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
analysis_name = "Interactions"

# Create directory if not exist
output_path = file.path(base_dir, comparison, analysis_name)
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Load the data
adata <- read_h5ad(setup$data_path)

# Define classification task
data <- adata[adata$obs[[setup$classification$output_column]]
              %in% setup$classification$comparison]

# Set intercept to EUR and any of the cancers
ancestry_column <- as.name(setup$classification$ancestry_column)
output_column <- as.name(setup$classification$output_column)
train_ancestry <- setup$classification$train_ancestry

# Relevel ancestry column to train ancestry as the intercept
data$obs <- data$obs |>
  mutate(
    {{ ancestry_column }} := factor(
      .data[[as.character(ancestry_column)]], # Access the column by its name
      levels = c(train_ancestry,
                 setdiff(unique(.data[[as.character(ancestry_column)]]),
                         train_ancestry))
    )
  )

# Find baseline
reference_levels <- sapply(data$obs, function(col) {
  if (is.factor(col)) levels(col)[1] else NA
})

baseline <- paste(reference_levels[[ancestry_column]],
                  reference_levels[[output_column]],
                  sep = ":")

# Assert the order of arrays
stopifnot(all(colnames(t(data$X)) ==
                row.names(data$obs[colnames(t(data$X)), ])))


# Filter genes 
# TODO - stronger filter maybe
keep <- filterByExpr(t(data$X), group=group)
data <- data[, keep]

# Create design matrix
# Dynamically construct the formula
design <- model.matrix(~ eval(ancestry_column) * eval(output_column),
                       data = data$obs)

# Remove "eval(ancestry_column)" and "eval(output_column)"
colnames(design) <- gsub("eval\\(ancestry_column\\)|eval\\(output_column\\)",
                         "",
                         colnames(design))

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
# TODO - Save the results from statistics

# Filter for interaction terms
limma_res_interactions <- filter(limma_res, str_detect(coef, ":"))

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
      x = coef,
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
      x = coef,
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









