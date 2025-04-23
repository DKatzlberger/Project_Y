# Remove start up messages
suppressPackageStartupMessages(
  {
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # DGE workflow and functional analysis
    library(edgeR)
    library(fgsea)
    # Parallelization
    library(parallel)
    # Visualization
    library(patchwork)
    library(ggrepel)
    library(ComplexHeatmap)
    library(GetoptLong)
    library(circlize)
    library(ggVennDiagram)
    # Standard libraries
    library(uuid)
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
    library(glue)
  }
)
# Custom functions
source("r_utils.R")
source("figure_themes.R")


# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if script is run with a command-line argument
if (length(args) > 0) {
  YAML_FILE <- args[1]
  # Check if it's a valid YAML file
  is_yaml_file(YAML_FILE)
} else {
  # Dev settings if no command-line argument provided
  cat("Running interactive mode for development. \n")
  YAML_FILE <- "example_settings_interactions.yaml"
  is_yaml_file(YAML_FILE)
}

# Input
# --- Settings
# Default
DEFAULT_FILE  <- "default_settings_interactions.yaml"
default_setup <- load_default_settings(DEFAULT_FILE)
# User
user_setup <- yaml.load_file(YAML_FILE)
# Default and user settings (user will overwrite default)
setup      <- modifyList(default_setup, user_setup)
# Check required settings
# Required settings
required_settings <- c(
  "output_column", "class_0", "class_1", "ancestry_column", 
  "train_ancestry", "infer_ancestry", "data_path", "tech", 
  "output_directory"
)
check_settings(setup, required_settings)
# Add info to settings
setup$date <- format(as.POSIXlt(Sys.time(), tz = "GMT"), "%Y-%m-%d %H:%M:%S") 
setup$id   <- toupper(substr(UUIDgenerate(), 1, 10))

# Create output directory
path_to_save_location <- setup$output_directory
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}
# Save the settings
save_name <- file.path(path_to_save_location, "Settings.yaml")
write_yaml(setup, save_name)

# --- Data
output_column     <- setup$output_column
class_0           <- setup$class_0
class_1           <- setup$class_1
ancestry_column   <- setup$ancestry_column
train_ancestry    <- setup$train_ancestry
infer_ancestry    <- setup$infer_ancestry

# Load data
data_path  <- setup$data_path
is_h5ad_file(data_path)
adata      <- read_h5ad(data_path)
# Check if columns exist in data
required_columns <- c(output_column, ancestry_column)
check_columns(adata$obs, required_columns)

# Define classification 
comparison <- c(class_0, class_1)
check_values(adata$obs, output_column, comparison)
adata      <- adata[adata$obs[[output_column]] %in% comparison]
# Factorize
adata$obs[[output_column]] <- factor(adata$obs[[output_column]], levels = comparison)

# Define ancestries
ancestries <- c(train_ancestry, infer_ancestry)
check_values(adata$obs, ancestry_column, ancestries)
adata      <- adata[adata$obs[[ancestry_column]] %in% ancestries]
# Factorize
adata$obs[[ancestry_column]] <- factor(adata$obs[[ancestry_column]], levels = ancestries)

# Interactions analysis
# Message
cat(sprintf("New analysis with id: %s; created: %s \n", setup$id, setup$date))
cat(sprintf("Save location: %s \n", path_to_save_location))
cat("-------------------------------------------------------------------- \n")

# Visualize: Sample sizes
save_name     <- file.path(path_to_save_location, "QC_sample_sizes.pdf")
# Plot
p_count       <- plot_output_column_count(adata$obs, ancestry_column, output_column)
p_proportions <- plot_output_column_proportion(adata$obs, ancestry_column, output_column)
# Combine
p <- p_count + p_proportions + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
save_ggplot(p, save_name, width = 6, height = 4)

# Add to settings
counts <- table(adata$obs[[output_column]], adata$obs[[ancestry_column]])
setup$n_class_0_train_ancestry <- counts[class_0, train_ancestry]
setup$n_class_1_train_ancestry <- counts[class_1, train_ancestry]
setup$n_class_0_infer_ancestry <- counts[class_0, infer_ancestry]
setup$n_class_1_infer_ancestry <- counts[class_1, infer_ancestry]


# --- Design matrix 
# Define groups to compare
adata$obs["group"] <- factor(
  paste(
    adata$obs[[ancestry_column]], 
    adata$obs[[output_column]], 
    sep = "."
  ), 
  levels = c(
    paste(train_ancestry, class_0, sep = "."),  
    paste(train_ancestry, class_1, sep = "."),  
    paste(infer_ancestry, class_0, sep = "."),    
    paste(infer_ancestry, class_1, sep = ".")     
  )
)
# Design formula and groups
formula <- as.formula(paste("~", "0 + group"))
groups  <- levels(adata$obs$group)
# Matrix
interaction_design <- model.matrix(~ 0 + group, data = adata$obs)
# Human and machine readable terms
colnames(interaction_design) <- gsub("group", "", colnames(interaction_design))
colnames(interaction_design) <- gsub("-", "_", colnames(interaction_design))

# --- Filter features 
# Settings
tech            <- setup$tech
data_type       <- setup$data_type
filter_features <- setup$filter_features

# Strength of filter
percentile <- setup$percentile

# Name of QC plot
save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
if (filter_features & tech == "transcriptomics"){

  # Print statement
  cat(sprintf("📊 Filter features for technology: %s. \n", tech))
  cat(sprintf("By count: counts -> norm factors -> logCPM -> count threshold (%s percentile) -> filter \n", percentile))

  # Transform to logCPM
  norm_factors <- calculate_tmm_norm_factors(adata$X)
  cpm_data     <- cpm(adata$X, norm_factors = norm_factors, log = TRUE)

  # Filter by signal/count
  min_counts        <- signal_by_percentile(cpm_data, percentile)
  filtered_features <- filter_by_signal(cpm_data, min_counts)
  # Subset
  filtered_data <- adata[, filtered_features]

  # Print statement
  cat(sprintf("By variance: counts -> norm factors -> logCPM -> variance threshold (%s percentile) -> filter \n", percentile))

  # Transform to logCPM
  norm_factors <- calculate_tmm_norm_factors(filtered_data$X)
  cpm_data     <- cpm(filtered_data$X, norm_factors = norm_factors, log = TRUE)

  # Filter by variance
  min_variance      <- variance_by_percentile(cpm_data, percentile)
  filtered_features <- filter_by_variance(cpm_data, var_threshold = min_variance)
  # Subset
  filtered_data <- filtered_data[, filtered_features]

  # Visualize: Filtering
  data_before <- log2(adata$X + 0.5)
  data_after  <- log2(filtered_data$X + 0.5)

  # Axis
  x_axis <- paste0("log2(", data_type, " + 0.5)")
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis)
  p_after  <- plot_mean_variance_trend(data_after, x_axis)
  # Combine
  p <- p_before / p_after
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (filter_features & tech == "methylation") {

  # Transform to mvalues
  mvals <- beta_to_mvalue(adata$X)
  
  # Filter by variance
  min_variance       <- variance_by_percentile(mvals, percentile)
  filtered_features  <- filter_by_variance(mvals, var_threshold = min_variance)
  # Subset
  filtered_data = adata[, filtered_features]
  
  # Visualize: Filtering
  data_before <- adata$X 
  data_after  <- filtered_data$X

  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis) 
  p_after  <- plot_mean_variance_trend(data_after, x_axis)
  # Combine
  p <- p_before / p_after
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (filter_features & tech == "proteomics"){

  # No filtering
  filtered_data = adata

  # Visualize: Filtering
  data_before <- adata$X 

  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis) 
  # Save
  save_ggplot(p_before, save_name, width = 6, height = 4)

} else{

  # No filtering
  filtered_data <- adata

  # Visualize: Filtering
  data_before <- adata$X 
  
  # Axis
  x_axis <- data_type
  # Plot
  p_before <- plot_mean_variance_trend(data_before, x_axis)
  # Save
  save_ggplot(p_before, save_name, width = 6, height = 4)

}
# Save feature 
save_name <- file.path(path_to_save_location, "Features.yaml")
write_yaml(filtered_data$var_names, save_name)
# Number of features
setup$n_features  <- ncol(filtered_data)
# Message
if (filter_features){
  cat(sprintf("Number of features after filtering: %s (%s). \n", ncol(filtered_data), ncol(adata)))
  cat("Check plot: 'QC_mean_variance_trend.pdf' for visualization. \n")
  cat("-------------------------------------------------------------------- \n")
}

# ---- Normalization/Transformation
# Settings
tech                 <- setup$tech
normalization        <- setup$normalization
normalization_method <- normalization_methods[[tech]][[normalization]]$"function"
values_output_name   <- normalization_methods[[tech]][[normalization]]$"output_name"
# Print statement
cat(sprintf("📊 Normalization of features for technology: %s. \n", tech))

# Transpose (rows = Genes, cols = Samples)
data_t <- t(filtered_data$X)

data_norm <- normalization_method(data_t, interaction_design)

# Extract normalized matrix (used for plotting)
data_norm_matrix <- if (is.list(data_norm) && !is.null(data_norm$E)) data_norm$E else data_norm

# Visualize: Normalization
# Density per sample
p_before <- plot_density_of_samples(filtered_data$X, x_axis_label = setup$data_type) 
p_after  <- plot_density_of_samples(t(data_norm_matrix), x_axis_label = values_output_name)
# Combine
p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Save
p_name <- file.path(path_to_save_location, "QC_density_normalized_values.pdf")
save_ggplot(p, p_name, width = 3, height = 3)

# Q-Q plots per gene
p_before <- plot_qq_of_genes(filtered_data$X, n_features = 5)
p_after  <- plot_qq_of_genes(t(data_norm_matrix), n_features = 5)
# Combine
p <- p_before / p_after 
# Save
save_name <- file.path(path_to_save_location, "QC_qq_normalized_values.pdf")
save_ggplot(p, save_name, width = 6, height = 4)

# Print statement 
cat("Check plot: 'QC_density_normalized_values.pdf' and 'QC_qq_normalized_values.pdf' for visualization. \n")
cat("-------------------------------------------------------------------- \n")

# --- Means model
# Print statement 
cat("📊 Means model summary: \n")
cat(sprintf("Formula: %s\n", deparse(formula)))
cat("Groups:\n")
cat(sprintf("  %s\n", paste(groups, collapse = paste0("\n  "))  ))

# Fit the model (means model)
limma_fit      <- lmFit(data_norm, design = interaction_design)
limma_fit      <- eBayes(limma_fit)
mean_model_res <- extract_results(limma_fit)

# Save
save_name <- file.path(path_to_save_location, "Limma_means.csv")
fwrite(mean_model_res, save_name)
cat("Check results: 'Limma_means.csv'. \n")
cat("-------------------------------------------------------------------- \n")

# --- Contrast 
# Terms
contrast_terms <- list(
  baseline_1      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}.Baseline"),
  baseline_2      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[2]}.Baseline"),
  relationship_1  = glue("{train_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relationship"),
  relationship_2  = glue("{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relatioship"),
  interaction     = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Interaction")
)

# Calculations
cols <- colnames(interaction_design)
contrast_calculations <- list(
  baseline_1      = glue("{cols[1]} - {cols[3]}"),
  baseline_2      = glue("{cols[2]} - {cols[4]}"),
  relationship_1  = glue("{cols[1]} - {cols[2]}"),
  relationship_2  = glue("{cols[3]} - {cols[4]}"),
  interaction     = glue("({cols[1]} - {cols[2]}) - ({cols[3]} - {cols[4]})")
)

# Print statement
cat("📊 Contrast summary: \n")
cat("Calculations: \n")
for (i in seq_along(contrast_calculations)) {
  cat(sprintf("  %-15s %s \n", names(contrast_calculations)[i], contrast_calculations[[i]]))
}

# Create contrast matrix
contrast_matrix <- makeContrasts(
  contrasts = contrast_calculations,
  levels    = interaction_design
)
colnames(contrast_matrix) <- contrast_terms

# Fit contrast
limma_fit_contrast <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit_contrast <- eBayes(limma_fit_contrast)
contrast_res       <- extract_results(limma_fit_contrast)

# Save
save_name <- file.path(path_to_save_location, "Limma_contrast.csv")
fwrite(contrast_res, save_name)
cat("Check results: 'Limma_contrast.csv'. \n")
cat("-------------------------------------------------------------------- \n")


# Save the settings
save_name <- file.path(path_to_save_location, "Settings.yaml")
write_yaml(setup, save_name)


# --- Visualize: Results
if (setup$visual_val){

  # Message
  cat("Visualizing results. \n")
  cat("-------------------------------------------------------------------- \n")
  # Create output directory
  path_to_save_location <- file.path(setup$output_directory, "Visual_val")
  if (!dir.exists(path_to_save_location)) {
    dir.create(path_to_save_location, recursive = TRUE)
  }

  # --- Baseline
  new_c <- c(ancestry_column, output_column, "comp_level")
  baseline_terms <- c(contrast_terms$baseline_1, contrast_terms$baseline_2)
  baseline       <- filter(contrast_res, coef %in% baseline_terms)
  baseline       <- separate(baseline, coef, into = new_c, sep = "\\.", remove = FALSE)
  # Add information
  baseline[[ancestry_column]] <- str_replace_all(baseline[[ancestry_column]], "_", " ")

  # Visualize: Baseline
  # Volcano plot
  p <- volcano_plot(
    baseline, 
    logFC_thr  = setup$logFC_thr,
    facet_cols = c(ancestry_column, output_column),
    point_size = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Baseline_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # MA plot
  p  <- ma_plot(
    baseline,
    logFC_thr    = setup$logFC_thr,
    y_axis_label = values_output_name,
    facet_cols   = c(ancestry_column, output_column),
    point_size   = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Baseline_ma.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # --- Relationship
  new_c <- c(ancestry_column, output_column, "comp_level")
  relationship_terms <- c(contrast_terms$relationship_1, contrast_terms$relationship_2)
  relationship       <- filter(contrast_res, coef %in% relationship_terms)
  relationship       <- separate(relationship, coef, into = new_c, sep = "\\.", remove = FALSE)
  # Add information
  relationship[[output_column]] <- str_replace_all(relationship[[output_column]], "_", " ")

  # Visualize: Relationship
  # Volcano plot
  p <- volcano_plot(
    relationship, 
    logFC_thr  = setup$logFC_thr,
    facet_cols = c(ancestry_column, output_column),
    point_size = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Relationship_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # MA plot
  p  <- ma_plot(
    relationship,
    logFC_thr    = setup$logFC_thr,
    y_axis_label = values_output_name,
    facet_cols   = c(ancestry_column, output_column),
    point_size   = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Relationship_ma.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)
  
  # --- Interaction
  new_c <- c(ancestry_column, output_column, "comp_level")
  interaction_term <- contrast_terms$interaction
  interaction      <- filter(contrast_res, coef %in% interaction_term)
  interaction      <- separate(interaction, coef, into = new_c, sep = "\\.", remove = FALSE)
  # Add information
  interaction[[ancestry_column]] <- str_replace_all(interaction[[ancestry_column]], "_", " ")
  interaction[[output_column]]   <- str_replace_all(interaction[[output_column]], "_", " ")

  # Visualize: Interaction
  # Volcano plot
  p <- volcano_plot(
    interaction, 
    logFC_thr  = setup$logFC_thr,
    facet_cols = c(ancestry_column, output_column),
    point_size = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Interaction_volcano.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # MA plot
  p  <- ma_plot(
    interaction,
    logFC_thr    = setup$logFC_thr,
    y_axis_label = values_output_name,
    facet_cols   = c(ancestry_column, output_column),
    point_size   = 1
  )
  # Save
  save_name <- file.path(path_to_save_location, "Interaction_ma.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

  # Significant interaction
  sig_interaction <- filter(interaction, adj.P.Val < 0.05 & abs(logFC) > setup$logFC_thr)
  sig_features    <- sig_interaction$Feature

  # Visualize: Sig interaction
  if (length(sig_features) != 0){
    
    # Number of features to visualize
    n_sig_features <- 10
    exp_list       <- list()
    for(feature in sig_features[1:n_sig_features]){
      # Meta data for each feature
      exp_list[[feature]] <- adata$obs |>
        mutate(zscore = scale(data_norm_matrix[feature,])) |>
        rownames_to_column("patient_id") |>
        remove_rownames()
    }
    # Features with normalised expression
    exp_zscore <- bind_rows(exp_list, .id = "Feature")

    # Boxplot
    p <- interactions_boxplot(
      exp_zscore, 
      x          = ancestry_column, 
      fill       = output_column,
      point_size = 0.8
      )
    # Save
    save_name <- file.path(path_to_save_location, "Interaction_boxplot.pdf")
    save_ggplot(p, save_name, width = 6, height = 4)

    # Heatmap
    interaction_heatmap(
      exp_zscore,
      output_column,
      ancestry_column,
      path_to_save_location
    )
  }
  cat(sprintf("Analysis %s finished.", setup$id))

} else{
  cat(sprintf("Analysis %s finished.", setup$id))
}






# # Results
# # Filter coeficients
# new_cols <- c("Ancestry", "Condition", "Difference")
# # Baseline
# baseline_terms <- c(contrast_terms$baseline_1, contrast_terms$baseline_2)
# baseline <- filter(contrast_res, coef %in% baseline_terms)
# baseline <- separate(baseline, coef, into = new_cols, sep = "\\.", remove = FALSE) |>
#   mutate(
#     Ancestry = toupper(Ancestry),
#     Condition = str_replace_all(Condition, "_", ""),
#     data_type = setup$data_type
#   ) 
# # Relationship
# relationship_terms <- c(contrast_terms$relationship_1, contrast_terms$relationship_2)
# relationship <- filter(contrast_res, coef %in% relationship_terms)
# relationship <- separate(relationship, coef, into = new_cols, sep = "\\.", remove = FALSE) |>
#   mutate(
#     Ancestry = toupper(Ancestry),
#     Condition = str_replace_all(Condition, "_", " "),
#     data_type = setup$data_type
#   ) 
# # Interaction
# interaction_term <- contrast_terms$interaction
# interaction <- filter(contrast_res, coef %in% interaction_term)
# interaction <- separate(interaction, coef, into = new_cols, sep = "\\.", fill = "right", remove = FALSE) |>
#   mutate(
#     Ancestry = toupper(Ancestry),
#     Condition = NA_character_,
#     Difference = "Interaction",
#     data_type = setup$data_type
#   ) 

# # Save 
# fwrite(baseline, file.path(path_to_save_location, "Baseline.csv"))
# fwrite(relationship, file.path(path_to_save_location, "Relationship.csv"))
# fwrite(interaction, file.path(path_to_save_location, "Interaction.csv"))

# # Results with significants
# logFC_thr <- 1
# sig_baseline <- filter(baseline, adj.P.Val < 0.05 & abs(logFC) > logFC_thr)
# sig_relationship <- filter(relationship, adj.P.Val < 0.05 & abs(logFC) > logFC_thr)
# sig_interaction <- filter(interaction, adj.P.Val < 0.05 & abs(logFC) > logFC_thr)

# # Results with trend
# trend_baseline <- filter(baseline, abs(logFC) > logFC_thr) |> 
#   group_by(Condition) |>
#   arrange(adj.P.Val, desc(abs(logFC))) |>
#   slice_head(n = 10)
# trend_relationship <- filter(relationship, abs(logFC) > logFC_thr) |> 
#   group_by(Ancestry) |>
#   arrange(adj.P.Val, desc(abs(logFC))) |>
#   slice_head(n = 10)
# trend_interactions <- filter(interaction, abs(logFC) > logFC_thr) |>
#   arrange(adj.P.Val, desc(abs(logFC))) |>
#   slice_head(n = 10)


# # Visualize
# print("Visualize.")
# # ---- Volcano plots ---- 
# # Baseline
# p <- volcano_plot(
#   baseline, 
#   logFC_thr, 
#   facet_rows = c("Difference", "Condition"),
#   facet_cols = "Ancestry"
#   )
# # Save
# save_name <- file.path(path_to_save_location, "Baseline_volcano.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)

# # Relationship
# p <- volcano_plot(
#   relationship, 
#   logFC_thr, 
#   facet_rows = c("Difference", "Condition"),
#   facet_cols = "Ancestry"
#   )
# # Save
# save_name <- file.path(path_to_save_location, "Relationship_volcano.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)

# # Interaction
# p <- volcano_plot(
#   interaction, 
#   logFC_thr, 
#   facet_rows = c("Difference"),
#   facet_cols = "Ancestry"
#   )
# # Save
# save_name <- file.path(path_to_save_location, "Interaction_volcano.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)

# # ---- Venn diagram ---- 
# venn_names <- c(comparison[1], comparison[2], "Interaction")
# venn_list <- setNames(vector("list", length(venn_names)), venn_names)

# # Baseline
# venn_baseline <- group_by(sig_baseline, Condition) |> 
#   summarise(Features = list(Feature)) |>
#   deframe()
# # Add to list
# venn_list[names(venn_baseline)] <- venn_baseline

# # Interaction
# venn_interaction <- group_by(sig_interaction, Difference)|> 
#   summarise(Features = list(Feature)) |>
#   deframe()
# # Add to list
# venn_list[names(venn_interaction)] <- venn_interaction

# # Replace emtpy lists
# venn_list <- lapply(venn_list, function(x) if (is.null(x)) character(0) else x)

# # Transform to venn readable data
# venn_data <- process_data(Venn(venn_list))
# # Plot
# p <- venn_diagram(venn_data)
# # Save
# save_name <- file.path(path_to_save_location, "Venn_diagram.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)

# # ---- Correlation relationship ----
# sig_baseline_features_1 <- filter(sig_baseline, Condition == comparison[1]) |> pull(Feature)
# sig_baseline_features_2 <- filter(sig_baseline, Condition == comparison[2]) |> pull(Feature)
# sig_interaction_features <- sig_interaction |> pull(Feature)

# # Sets
# set_names <- c(
#   paste("Baseline:", comparison[1]),
#   paste("Baseline:", comparison[2]),
#   paste("Interaction + Baseline:", comparison[1]),
#   paste("Interaction + Baseline:", comparison[2]),
#   "Interaction",
#   "Non-significant"
# )

# # Colors an fill
# set_colors <- c(
#   "black",  
#   "black", 
#   "black", 
#   "black", 
#   "black", 
#   "grey"
# )

# set_fill <- c(
#   "#a3ddcb",  
#   "#e7c6fb",  
#   "#027c58",  
#   "#9c06f9",  
#   "red",  
#   "grey"
# )

# # Rename
# names(set_colors) <- set_names
# names(set_fill) <- set_names

# # Prepare data
# p_data <- relationship |> 
#   select(-where(is.numeric), -coef, logFC) |>
#   pivot_wider(
#     names_from = Ancestry,
#     values_from = logFC
#   ) |>
#   mutate(
#     with_baseline_1   = Feature %in% sig_baseline_features_1,
#     with_baseline_2   = Feature %in% sig_baseline_features_2,
#     with_interaction  = Feature %in% sig_interaction_features, 
#     coloring          = case_when(
#       !with_interaction & with_baseline_1  ~ set_names[1], 
#       !with_interaction & with_baseline_2  ~ set_names[2],  
#       with_interaction  & with_baseline_1  ~ set_names[3],
#       with_interaction  & with_baseline_2  ~ set_names[4],
#       with_interaction  & !with_baseline_1 ~ set_names[5],    
#       with_interaction  & !with_baseline_2 ~ set_names[5],                
#       TRUE                                 ~ set_names[6]
#       )
#   )

# # Threshold line
# lm_formula <- as.formula(paste(toupper(inf_ancestry), "~", toupper(train_ancestry)))
# lm_model <- lm(lm_formula, data = p_data)
# # Add threshold
# residual_thr <- 1
# upper_threshold <- fitted(lm_model) + 1
# lower_threshold <- fitted(lm_model) - 1

# # Correlation 
# cor_pearson <- cor(p_data[toupper(train_ancestry)], p_data[toupper(inf_ancestry)], method = "pearson")
# cor_spearman <- cor(p_data[toupper(train_ancestry)], p_data[toupper(inf_ancestry)], method = "spearman")
# # Label
# pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
# spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

# # Axis name
# x_axis <- paste(toupper(train_ancestry), unique(p_data$Condition), "logFC")
# y_axis <- paste(toupper(inf_ancestry), unique(p_data$Condition), "logFC")

# # Plot
# p <- p_data |>
#   ggplot(
#     aes(
#       x = !!sym(toupper(train_ancestry)),
#       y = !!sym(toupper(inf_ancestry))
#     )
#   ) +
#   geom_point(
#     data = p_data |> filter(coloring == "Non-significant"),
#     aes(
#       fill = coloring,
#       color = coloring
#       ),
#     shape = 21,
#     size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     data = p_data |> filter(coloring %in% c(
#       set_names[1],
#       set_names[2]
#       )
#     ),
#     aes(
#       fill = coloring,
#       color = coloring
#       ),
#     shape = 21,
#     size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     data = p_data |> filter(coloring %in% c(
#       set_names[3],
#       set_names[4]
#       )
#     ),
#     aes(
#       fill = coloring,
#       color = coloring
#       ),
#     shape = 21,
#     size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     data = p_data |> filter(coloring %in% c(
#       set_names[5]
#       )
#     ),
#     aes(
#       fill = coloring,
#       color = coloring
#       ),
#     shape = 21,
#     size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     aes(
#       fill = "dummy1",
#       color = "dummy1"
#       ),
#       shape = 21,
#       size = 0,
#       stroke = 0,
#       alpha = 0
#   ) +
#   geom_point(
#     aes(
#       fill = "dummy2",
#       color = "dummy2"
#       ),
#       shape = 21,
#       size = 0,
#       stroke = 0,
#       alpha = 0
#   ) +
#   geom_smooth(
#     method = "lm",
#     color = "blue",
#     linewidth = (0.5 / 2),
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(y = upper_threshold),
#     color = "blue", 
#     linetype = "dashed",
#     linewidth = (0.5 / 1.5),
#     alpha = 0.5
#   ) + 
#     geom_line(
#     aes(y = lower_threshold), 
#     color = "blue", 
#     linetype = "dashed",
#     linewidth = (0.5 / 1.5),
#     alpha = 0.5
#   ) + 
#   scale_fill_manual(
#     name = "",
#     values = c(set_fill, "dummy1" = "white", "dummy2" = "white"),
#     breaks = c(set_names, "dummy1", "dummy2"),
#     labels = c(set_names, pearson_label, spearman_label)
#   ) +
#   scale_color_manual(
#     name = "",
#     values = c(set_colors, "dummy1" = "white", "dummy2" = "white"),
#     breaks = c(set_names, "dummy1", "dummy2"),
#     labels = c(set_names, pearson_label, spearman_label)
#   ) +
#   labs(
#     x = x_axis,
#     y = y_axis
#   ) + 
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_white_strip() +
#   theme_small_legend() +
#   theme(
#     legend.position = c(0.05, 0.95),       
#     legend.justification = c(0, 1),  
#     legend.background = element_rect(fill = "transparent", color = NA)
#   )

# # Save
# save_name <- file.path(path_to_save_location, "Relationship_correlation.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)

# # ---- Quality of DGE -----
# deg_interaction <- interaction |>
#   mutate(
#     Significance = case_when(
#       adj.P.Val < 0.05 & abs(logFC) > logFC_thr & logFC > 0 ~ "Up-regulated",
#       adj.P.Val < 0.05 & abs(logFC) > logFC_thr & logFC < 0 ~ "Down-regulated",
#       TRUE                                                  ~ "Non-significant"
#     )
#   )

# # Set names
# set_names <- c(
#   "Up-regulated",
#   "Down-regulated",
#   "Non-significant"
# )

# # Colors
# set_colors <- c(
#   "black",  
#   "black", 
#   "grey"
# )

# # Fill
# set_fill <- c(
#   "#e7e705",
#   "#02b102",
#   "grey"
# )

# # Rename
# names(set_colors) <- set_names
# names(set_fill) <- set_names

# p <- deg_interaction |>
#   ggplot(
#     aes(
#       x = AveExpr,
#       y = logFC
#     )
#   ) +
#   geom_point(
#     data = filter(deg_interaction, Significance == "Non-significant"),
#     aes(
#       fill = Significance,
#       color = Significance
#     ),
#     shape = 21,
#     size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     data = filter(deg_interaction, Significance == "Up-regulated"),
#     aes(
#       fill = Significance,
#       color = Significance,
#       size = -log(adj.P.Val)
#     ),
#     shape = 21,
#     #size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_point(
#     data = filter(deg_interaction, Significance == "Down-regulated"),
#     aes(
#       fill = Significance,
#       color = Significance,
#       size = -log(adj.P.Val)
#     ),
#     shape = 21,
#     #size = 1.5,
#     stroke = 0.1
#   ) +
#   geom_hline(
#     yintercept = 0, 
#     color = "black", 
#     linewidth = (0.5 / 2)
#   ) +
#   scale_fill_manual(
#     name = "Interaction",
#     values = c(set_fill),
#     breaks = set_names
#   ) +
#   scale_color_manual(
#     name = "Interaction",
#     values = c(set_colors),
#     breaks = set_names
#   ) +
#   scale_size_binned(
#     range = c(0.5, 3)
#   ) +
#   labs(
#     x = paste("Average", values_output_name)
#   ) +
#   guides(
#     size = guide_legend(
#       title.position = "top",  
#       direction = "horizontal"
#     )
#   ) +
#   theme_nature_fonts() +
#   theme_white_background () +
#   theme_small_legend() +
#   theme(
#     legend.position = c(0.95, 0.05),       
#     legend.justification = c("right", "bottom"),  
#     legend.background = element_rect(fill = "transparent", color = NA)
#   )

# # Save
# save_name <- file.path(path_to_save_location, "Interaction_MA.pdf")
# save_ggplot(p, save_name, width = 3, height = 3)


# if (nrow(sig_interaction) > 0){
# # ---- Interaction Heatmap ----
# exp_list <- list()
# for(gg in sig_interaction$Feature){
#   exp_list[[gg]] <- adata$obs |>
#     mutate(z_score = scale(data_norm_matrix[gg,])) |>
#     rownames_to_column("patient_id") |>
#     remove_rownames()
# }
# # Expression zscore
# expression_zscore <- bind_rows(exp_list, .id = "Feature")

# # Expression matrix
# expression_matrix <- expression_zscore |>
#   select(Feature, patient_id, z_score) |>
#   spread(key = patient_id, value = z_score) |>
#   column_to_rownames("Feature") |>
#   as.matrix()

# # Condition colors
# condition_colors <- setNames(c("#027c58", "purple"), comparison)
# genetic_ancestry_colors <- c(
#   "admix"= "#ff4d4d", 
#   "afr" = "#ff9900", 
#   "amr" =  "#33cc33",
#   "eur" = "#3399ff", 
#   "eas" = "#cc33ff", 
#   "sas" = "#ffcc00"
#   )

# # Annotation
# annotations <- expression_zscore |>
#   select(patient_id, all_of(c(output_column, ancestry_column))) |>
#   distinct() |>
#   arrange(match(patient_id, colnames(expression_zscore)))

# annotations <- annotations |>
#   mutate(
#     group = paste(annotations[[ancestry_column]], annotations[[output_column]], sep = "."),
#     group = factor(group, levels = str_replace(colnames(interaction_design), "_", "-")),
#   ) 

# # Heatmap annotation
# heatmap_annotation <- HeatmapAnnotation(
#   # Annotations
#   Ancestry = anno_empty(
#     border = FALSE, 
#     height = unit(0.2, "cm"), 
#     show_name = TRUE
#   ),
#   Condition = annotations[[output_column]], 
#   # Colors
#   col = list(
#     Condition = condition_colors
#   ), 
#   # Fonts
#   annotation_name_gp = gpar(fontsize = 5),
#   border = FALSE,
#   show_annotation_name = TRUE,
#   annotation_name_side = "left",
#   # Height
#   simple_anno_size = unit(0.2, "cm"),
#   # Legend
#   show_legend = FALSE
# )

# # Legends
# ancetsry_lgd <- Legend(
#   title = "Ancestry", 
#   at = c(
#     toupper(train_ancestry), 
#     toupper(inf_ancestry)
#     ), 
#   legend_gp = gpar(
#     fill = c(
#       genetic_ancestry_colors[train_ancestry], 
#       genetic_ancestry_colors[inf_ancestry]
#       )
#     ),
#   # Font
#   labels_gp = gpar(fontsize = 5),
#   title_gp = gpar(fontsize = 5, fontface = "plain"),
#   # Size
#   grid_height = unit(0.3, "cm"), 
#   grid_width = unit(0.3, "cm"),
#   # Direction
#   direction = "horizontal"
#   )

# condition_lgd <- Legend(
#   title = "Condition", 
#   at = c(
#     comparison[1], 
#     comparison[2]
#     ), 
#   legend_gp = gpar(
#     fill = c(
#       condition_colors[comparison[1]], 
#       condition_colors[comparison[2]]
#       )
#     ),
#   # Font
#   labels_gp = gpar(fontsize = 5),
#   title_gp = gpar(fontsize = 5, fontface = "plain"),
#   # Size
#   grid_height = unit(0.3, "cm"), 
#   grid_width = unit(0.3, "cm"),
#   # Direction
#   direction = "horizontal"
#   )

# lgd_list <- c(ancetsry_lgd, condition_lgd)

# # Heatmap
# heatmap <- Heatmap(
#   expression_matrix,
#   name = "z-score", 
#   show_column_names = FALSE,
#   # Clustering
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   # Split
#   column_split = annotations$group,
#   cluster_column_slices = FALSE,
#   column_title = NULL,
#   # Annotations
#   top_annotation = heatmap_annotation,
#   # Fonts
#   row_names_gp = gpar(fontsize = 5),
#   row_title_gp = gpar(fontsize = 5),
#   row_names_side = "left",
#   # Legends
#   heatmap_legend_param = list(
#     title_gp = gpar(fontsize = 5, fontface = "plain"), 
#     labels_gp = gpar(fontsize = 5, fontface = "plain"),
#     grid_height = unit(0.3, "cm"), 
#     grid_width = unit(0.3, "cm"),
#     direction = "horizontal"
#     )
#   )

# # Save
# save_name <- file.path(path_to_save_location, "Interaction_heatmap.pdf")
# pdf(save_name, width = 6 , height = 3)

# # Add legend
# draw(
#   heatmap, 
#   # Legend
#   annotation_legend_list = lgd_list,
#   heatmap_legend_side = "bottom",
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
#   )

# # Spanning groups
# group_block_anno = function(group, empty_anno, gp = gpar(), 
#     label = NULL, label_gp = gpar()) {

#     seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
#     loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
#     seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
#     loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

#     seekViewport("global")
#     grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
#         just = c("left", "bottom"), gp = gp)
#     if(!is.null(label)) {
#         grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
#     }
# }

# # Add the blocks
# group_block_anno(
#   1:2, "Ancestry", 
#   gp = gpar(
#     fill = genetic_ancestry_colors[train_ancestry],
#     col = NA
#     )
#   )
# group_block_anno(
#   3:4, "Ancestry", 
#   gp = gpar(
#     fill = genetic_ancestry_colors[inf_ancestry],
#     col = NA
#     )
#   )

# # Close pdf and delete heatmap object
# dev.off()  
# rm(heatmap, envir = .GlobalEnv)

# # ---- Interaction Boxplots -----
# exp_list <- list()
# for(gg in head(sig_interaction$Feature, 10)){
#   exp_list[[gg]] <- adata$obs |>
#     mutate(z_score = scale(data_norm_matrix[gg,])) |>
#     rownames_to_column("patient_id") |>
#     remove_rownames()
# }
# # Expression zscore
# expression_zscore <- bind_rows(exp_list, .id = "Feature")
# # Number or rows
# n_facets <- length(unique(expression_zscore$Feature))
# n_col = 5
# n_row <- floor(n_facets / n_col)  
#   if (n_facets %% n_col != 0) {    
#     n_row <- n_row + 1
#   }

# p <- interactions_boxplot(
#   expression_zscore, 
#   ancestry_column, 
#   output_column,
#   nrow = n_row,
#   ncol = n_col
#   )
# # Save
# save_name <- file.path(path_to_save_location, "Interaction_boxplots.pdf")
# save_ggplot(p, save_name, width = 6, height = 3)


# } else{
#   print("Maybe visualize trends.")
# }



# # Functional analysis
# database <- "data/downloads/geneset-libraries/MSigDB_Hallmark_2020.txt"
# database <- read_enrichR_database(database)

# # Baseline
# baseline_enrich <- baseline |> 
#   group_split(Condition) |> 
#   map(
#     function(x){
#       enrichment <- perform_gsea(x, database, rank_by = "logFC") |>
#       mutate(
#         Ancestry    = unique(x$Ancestry),
#         Condition   = unique(x$Condition),
#         Difference  = unique(x$Difference),
#         data_type   = setup$data_type
#         )
#     }
#   ) |>
#   bind_rows()

# # Relationship
# relationship_enrich <- relationship |> 
#   group_split(Ancestry) |> 
#   map(
#     function(x){
#       enrichment <- perform_gsea(x, database, rank_by = "logFC") |>
#       mutate(
#         Ancestry    = unique(x$Ancestry),
#         Condition   = unique(x$Condition),
#         Difference  = unique(x$Difference),
#         data_type   = setup$data_type
#         )
#     }
#   ) |>
#   bind_rows()

# # Interaction
# interaction_enrich <- perform_gsea(interaction, database, rank_by = "logFC") |>
#   mutate(
#     Ancestry    = unique(interaction$Ancestry),
#     Condition   = NA_character_,
#     Difference  = unique(interaction$Difference),
#     data_type   = setup$data_type
#   )

# # Save
# fwrite(baseline_enrich, file.path(path_to_save_location, "Baseline_enrichment.csv"))
# fwrite(relationship_enrich, file.path(path_to_save_location, "Relationship_enrichment.csv"))
# fwrite(interaction_enrich, file.path(path_to_save_location, "Interaction_enrichment.csv"))

# # Visualize
# # ---- Dot plot ----
# # Baseline
# save_name <- file.path(path_to_save_location, "Baseline_enrichement.pdf")
# p <- fgsea_plot(
#   baseline_enrich,
#   x = "Condition",
#   facet_cols = "Ancestry"
#   )
# # Save
# save_ggplot(p, save_name, width = 3, height = 5)

# # Relationships
# save_name <- file.path(path_to_save_location, "Relationship_enrichement.pdf")
# p <- fgsea_plot(
#   relationship_enrich,
#   x = "Ancestry",
#   facet_cols = "Condition"
#   )
# # Save
# save_ggplot(p, save_name, width = 3, height = 5)

# # Interactions
# save_name <- file.path(path_to_save_location, "Interaction_enrichement.pdf")
# p <- fgsea_plot(
#   interaction_enrich,
#   x = "Difference",
#   facet_cols = "Ancestry"
#   )
# # Save
# save_ggplot(p, save_name, width = 3, height = 5)


# # Cluster profiler 
# GO_annotations <- fread("data/downloads/geneset-libraries/GO_annotation.csv")
# GO_ids <- inner_join(interactions, GO_annotations, by=c("Feature"="gene_name"))
# # Filter significant genes
# GO_sig_interactions <- filter(GO_ids, adj.P.Val < 0.05) |> pull(gene_id)
# # Background (EUR condition1 vs condition2)
# GO_background <- GO_ids |> pull(gene_id)
# GO_background <- setdiff(GO_background, GO_sig_interactions)
# # Remove input
# GO_background <- setdiff(GO_background, GO_sig_interactions)

# # Overrepresentation
# GO_overrepresentation <- enrichGO(
#   gene = GO_sig_interactions,
#   universe = GO_background,
#   keyType = "ENSEMBL",
#   OrgDb = org.Hs.eg.db, 
#   ont = "ALL",
#   pAdjustMethod = "BH", 
#   qvalueCutoff = 0.05, 
#   readable = TRUE)
# GO_results  <- as_tibble(GO_overrepresentation) |>
#   mutate(Ancestry = toupper(inf_ancestry))
# # Save
# fwrite(GO_results, file.path(path_to_save_location, "Interactions_overrepresentation.csv"))



# ---- Correlation of baselines (baseline_scatter_plot) -------
# # Reformat data
# wide_baseline <- baseline |> 
#   mutate(Comparison = paste(Comparison, Condition)) |>
#   dplyr::select(Feature, Comparison, logFC) |>
#   pivot_wider(names_from = Comparison, values_from = logFC) |>
#   rename_with(~ str_replace_all(., " ", "_"))

# # Fit linear model and get predictions
# x_col <- names(wide_baseline)[2]
# y_col <- names(wide_baseline)[3]
# lm_model <- lm(as.formula(paste(y_col, "~", x_col)), data = wide_baseline)

# # Coefficients
# intercept <- coef(lm_model)[1]
# slope <- coef(lm_model)[2]

# predicted_values <- fitted(lm_model)

# residuals <- residuals(lm_model)
# abs_residuals <- abs(residuals)

# # Define outliers
# residual_threshold <- 1
# upper_threshold <- predicted_values + residual_threshold
# lower_threshold <- predicted_values - residual_threshold

# # Calculate the correlation
# cor_pearson <- cor(wide_baseline[[x_col]], wide_baseline[[y_col]], method = "pearson")
# cor_spearman <- cor(wide_baseline[[x_col]], wide_baseline[[y_col]], method = "spearman")
# # Label
# label_pearson <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
# label_spearman <- bquote(R[Spearman] == .(round(cor_pearson, 3)))

# # Reformat data
# scatter_plot_data <- baseline |>
#   mutate(Comparison = paste(Comparison, Condition, "(logFC)")) |>
#   dplyr::select(Feature, Comparison, logFC) |>
#   pivot_wider(names_from = Comparison, values_from = logFC) |>
#   mutate(
#     abs_residual = abs_residuals,
#     above_threshold = abs_residual > residual_threshold,
#     upper_threshold = upper_threshold,
#     lower_threshold = lower_threshold
#   ) |>
#   mutate(
#     with_interactions = Feature %in% sig_interactions, 
#     with_baseline = Feature %in% sig_baseline,
#     coloring = case_when(
#       with_interactions & with_baseline ~ "Interaction + Baseline",
#       with_interactions &! with_baseline ~ "Interaction",    
#       with_baseline & !with_interactions ~ "Baseline",
#       TRUE ~ "Rest"                                
#     )
#   )

# # Residuals above threshold
# sig_residuals <- scatter_plot_data |>
#   filter(abs(residuals) > logFC_threshold) |>
#   pull(Feature)

# # Legend 
# # Legend x-position
# max_name_length <- max(nchar(unique(scatter_plot_data$coloring)))
# legend_x_pos <- min(0.5, 0.05 + 0.005 * max_name_length)
# # Legend y-position
# max_items <- length(unique(scatter_plot_data$coloring))
# legend_y_pos <- max(0.3, 0.9 - 0.01 * max_items)

# # Plot
# baseline_scatter_plot <- scatter_plot_data %>%
#   ggplot(
#     aes(
#       x = .data[[names(.)[2]]],
#       y = .data[[names(.)[3]]],
#       color = coloring
#     )
#   ) +
#   geom_point(
#     data = scatter_plot_data |> filter(coloring == "Rest"),
#     aes(color = coloring),
#     size = point_size,
#     show.legend = FALSE
#   ) +
#   geom_point(
#     data = scatter_plot_data |> filter(coloring == "Rest"),
#     aes(color = "dummy"),
#     size = NA,   
#     show.legend = FALSE 
#   ) +
#   geom_point(
#     data = scatter_plot_data |> filter(coloring == "Baseline"),
#     aes(color = coloring),
#     size = point_size
#   ) +
#   geom_point(
#     data = scatter_plot_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
#     aes(color = coloring),
#     size = point_size
#   ) +
#   geom_smooth(
#     method = "lm", 
#     color = "blue",
#     linewidth = point_size,
#     alpha = 0.5
#   ) + 
#   geom_line(
#     aes(y = upper_threshold), 
#     linetype = "dashed", 
#     color = "blue",
#     linewidth = point_size,
#     alpha = 0.5
#   ) +  
#   geom_line(
#     aes(y = lower_threshold), 
#     linetype = "dashed", 
#     color = "blue",
#     linewidth = point_size,
#     alpha = 0.5
#   ) +
#   scale_color_manual(
#     values = c(
#       "Interaction + Baseline" = "red",
#       "Interaction" = "darkgreen", 
#       "Baseline" = "gold", 
#       "Rest" = "lightgrey",
#       "dummy" = "black"
#     ),
#     breaks = c(
#       "Interaction + Baseline", 
#       "Interaction", 
#       "Baseline",
#       "Rest",
#       "dummy"
#     ),
#     labels = c(
#       "Interaction + Baseline", 
#       "Interaction", 
#       "Baseline",
#       label_pearson,
#       label_spearman
#     )
#   ) +
#   scale_x_continuous(
#     limits = c(
#       min(scatter_plot_data[[names(scatter_plot_data)[2]]]) - 1, 
#       max(scatter_plot_data[[names(scatter_plot_data)[2]]]) + 1
#     )  
#   ) +
#   scale_y_continuous(
#     limits = c(
#       min(scatter_plot_data[[names(scatter_plot_data)[3]]]) - 1, 
#       max(scatter_plot_data[[names(scatter_plot_data)[3]]]) + 1
#     )  
#   ) +
#   geom_text_repel(
#     data = scatter_plot_data |> filter(Feature %in% top_interactions),
#     aes(label = Feature),
#     size = 1.5,
#     color = "black",  
#     segment.color = "black",  
#     min.segment.length = 0
#   ) +
#   theme_nature_fonts() +
#   theme_white_background() +
#   theme_small_legend() +
#   theme(
#     legend.position = c(legend_x_pos, legend_y_pos), 
#     legend.title = element_blank(),  
#     legend.spacing.x = unit(0, "cm"), 
#     legend.spacing.y = unit(0, "cm"),
#     legend.background = element_rect(fill = "transparent"),  
#     legend.key.spacing = unit(0, "cm")
#   )

# # Save the plot
# ggsave(filename = "Baseline_scatter_plot.pdf",
#        plot = baseline_scatter_plot,
#        path = path_to_save_location,
#        height = 3, width = 3)




# # Functional Analysis
# # GSEA
# # Rank the genes based on (i) logFC (ii) pvalue
# # Database to enrich
# databases <- c("data/downloads/geneset-libraries/MSigDB_Hallmark_2020.txt",
#                "data/downloads/geneset-libraries/Reactome_Pathways_2024.txt",
#                )

# # Store all database terms
# fgsea_enrichment <- data.frame()
# fgsea_plots <- list()
# for (database in databases){
#   # Run fgsea
#   fgsea_logFC <- perform_gsea(interactions, database, rank_by = "logFC") |>
#     arrange(desc(abs(NES))) |>
#     mutate(Top50 = row_number() <= 50,
#            Top30 = row_number() <= 30,
#            Top10 = row_number() <= 10)

#   # Database name
#   db_name_ <- sub("\\.txt$", "", basename(database)) 
#   db_name <- gsub("_", " ", db_name_)

#   # Combine results
#   combined_fgsea <- fgsea_logFC |>
#     mutate(Database = db_name)
  
#   # Add to the list
#   fgsea_enrichment <- bind_rows(fgsea_enrichment, combined_fgsea)
  
#   # Visualize
#   top <- combined_fgsea |>
#     filter(Top30 == TRUE)
  
#   # Generate the GSEA plot
#   fgsea_plot <- top |>
#     ggplot(
#       aes(
#         x = ranked_by,
#         y = pathway,
#         color = NES,
#         size = pmin(-log10(padj), 5)
#         )
#     ) +
#     geom_point() +
#     scale_size_binned(
#      range = c(1, 3)    
#     ) +
#     scale_color_gradient2(
#       high = "red", 
#       mid = "white", 
#       low = "blue"
#     ) +
#     labs(
#       title = db_name,
#     )
#   # Add the plot to the list
#   fgsea_plots[[db_name]] <- fgsea_plot

#   # Save the plot
#   ggsave(
#     filename = paste0("Interactions_FGSEA_", db_name_, ".pdf"),
#     plot = fgsea_plot,
#     path = path_to_save_location,
#     width = 10,                                          
#     height = 5
#   )
# }

# # Save fgsea dataframe
# fgsea_enrichment <- fgsea_enrichment |>
#   mutate(Ancestry = toupper(setup$classification$infer_ancestry),
#          Interaction = paste(toupper(train_ancestry), "vs", toupper(setup$classification$infer_ancestry)))
# fwrite(fgsea_enrichment, file.path(path_to_save_location, "Interactions_FGSEA_enrichment.csv"))


# # clusterProfiler
# # Background genes comparison between cancers in European -> will look at overrepresented genes compared to European 
# # Substract the input from the background
# GO_annotations <- fread("data/downloads/geneset-libraries/GO_annotation.csv")
# # Merge results with annotations (will lose genes that are not in GO terms)
# GO_ids <- inner_join(interactions, GO_annotations, by=c("Feature"="gene_name"))
# # Define backround and significant genes
# GO_all_genes <- GO_ids$gene_id
# GO_sig_interactions <- filter(GO_ids, adj.P.Val < 0.05)
# GO_sig_genes <- GO_sig_interactions |> pull(gene_id)
# # # GO enrichment 
# # Background genes comparison between cancers in European
# GO_enrichment <- enrichGO(gene = GO_sig_genes,
#                           universe = GO_all_genes,
#                           keyType = "ENSEMBL",
#                           OrgDb = org.Hs.eg.db, 
#                           ont = "ALL", 
#                           pAdjustMethod = "BH", 
#                           qvalueCutoff = 0.05, 
#                           readable = TRUE)

# GO_results  <- as_tibble(GO_enrichment)
# # dotplot(GO_enrichment, showCategory=50)

# # # Visualize:
# GO_plot <- GO_results |>
#   separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) |>
#   mutate(GeneRatio = numerator / denominator)  |>
#   dplyr::select(-numerator, -denominator) |>
#   ggplot(aes(
#     x = GeneRatio,
#     y = Description,
#     size = -log(p.adjust),
#     color = FoldEnrichment
#   )
#   ) +
#   geom_point() +
#   scale_color_gradient(low = "blue", high = "red") 

# # Save
# ggsave(
#     filename = "Interactions_GO_overrepresentation.pdf",
#     plot = GO_plot,
#     path = path_to_save_location,
#     width = 12,                                          
#     height = 8
#   )

# # Netplot
# sig_foldchanges <- sig$logFC
# names(sig_foldchanges) <- sig$Feature
# cnetplot(GO_enrichment, 
#          categorySize="pvalue", 
#          showCategory = 5, 
#          foldChange=sig_foldchanges, 
#          vertex.label.font=6)

# res_entrez <- filter(GO_ids, entrezid != "NA")
# res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
# foldchanges <- res_entrez$logFC
# names(foldchanges) <- res_entrez$entrezid
# foldchanges <- sort(foldchanges, decreasing = TRUE)
# gseaKEGG <- gseKEGG(geneList = foldchanges, 
#               organism = "hsa", 
#               minGSSize = 20, 
#               pvalueCutoff = 0.05, 
#               verbose = FALSE)


# Find ancestry-specific enricht terms 
# Databases for gsea
# db <- read_enrichR_database("data/geneset-libraries/MSigDB_Hallmark_2020.txt")
# # 1. Filter for interaction effects
# # 2. Split 'interactions' by ancestry
# # 3. Run FGSEA
# interactions <- limma_fit_res |> filter(str_detect(coef, ":"))
# ancestry_specific_interactions <- interactions |> group_split(coef) 

# # Function to run FGSEA
# # Rank based on logFC
# run_fgsea <- function(split_df, gene_sets) {
#   # Create ranked list of genes
#   ranked_genes <- split_df |>
#     arrange(desc(logFC)) |>  # Order by logFC
#     distinct(Feature, .keep_all = TRUE) |>  # Ensure unique genes
#     pull(logFC, Feature)  # Pull named vector: logFC values with gene names
  
#   # Perform FGSEA
#   fgsea_results <- fgsea(pathways = gene_sets,
#                          stats = ranked_genes
#                          ) 
  
#   # Return results as a tibble
#   return(as_tibble(fgsea_results))
# }

# # Apply FGSEA to each split and store results
# fgsea_results_list <- lapply(ancestry_specific_interactions, function(split) {
#   coef_name <- unique(split$coef) # Name of the group
#   results <- run_fgsea(split, db)
#   results <- results |> mutate(coef = coef_name) # Add group info
#   return(results)
# })

# # Bind list object into dataframe
# fgsea_results <- bind_rows(fgsea_results_list)

# # Visualize 
# fgsea_results |>
#   mutate(coef = toupper(str_remove(coef, ":.*"))) |>
#   ggplot(
#     aes(
#       x = NES,
#       y = pathway,
#       size = -log10(padj)
#     )
#   ) + 
#   geom_point() +
#   facet_grid(
#     cols = vars(coef)
#   )











# Ranked gene list based on t-value
# ranked_interactions <- setNames(ancestry_specific_interactiona[[1]]$t, ancestry_specific_interactiona[[1]]$Feature)
# # Results
# res <- fgsea(
#   pathways = db,
#   stats = ranked_interactions)



# # Open pdf file
# pdf(file.path(output_path, "Voom_plot.pdf"))
# par(mfrow=c(1,2))
# # Normalization (logCPM)
# data_norm <- voom(t(data$X), design, plot = TRUE)
# # Fit the model
# limma_fit <- lmFit(data_norm, design = design)
# # Ebayes
# limma_fit <- eBayes(limma_fit)
# # Plot
# plotSA(limma_fit, main="Final model: Mean-variance trend")
# # Close the PDF device
# dev.off()

# # Extract results
# limma_res <- list()
# for (x in colnames(coef(limma_fit))){
#   # Extract for each coefficient all genes
#   limma_res[[x]] <- topTable(limma_fit, coef = x, number = Inf) |>
#      rownames_to_column("Feature")
# }
# limma_res <- bind_rows(limma_res, .id = "coef")
# limma_res <- filter(limma_res, coef != "(Intercept)")

# # Save the results from statistics
# fwrite(limma_res, 
#        file.path(output_path, "Interaction_results.csv"))

# # Filter for interaction terms
# limma_res_interactions <- filter(limma_res, str_detect(coef, ":"))

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
# volcano_plot <- limma_res_interactions |>
#   ggplot(
#     aes(
#       x = logFC,
#       y = -log10(adj.P.Val),
#       color = (adj.P.Val < 0.05 & abs(logFC) > 2)
#     )
#   ) +
#   geom_point(alpha = 0.8, size = 3) +
#   facet_grid(cols = vars(coef)) + 
#   geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
#   scale_color_manual(
#     values = c("TRUE" = "red", "FALSE" = "grey"),
#     name = 'Above thresholds'
#   ) 

# # Save the volcano plot
# ggsave(file.path(output_path, "Volcano_plot.pdf"),
#        plot = volcano_plot, height = 6, width = 12)

# # Pvalue distribution
# p_value_dist <- limma_res |>
#   ggplot(aes(x = P.Value,
#   fill=factor(floor(AveExpr)))) +
#   geom_histogram() +
#   facet_grid(rows = vars(coef))

# # Save the distribution plot
# ggsave(file.path(output_path, 
#        "P_value_plot.pdf"),
#        plot = p_value_dist, height = 10, width = 6)

# # Plot for interactions
# p_value_dist <- limma_res_interactions |>
#   ggplot(aes(x = P.Value,
#   fill=factor(floor(AveExpr)))) +
#   geom_histogram() +
#   facet_grid(rows = vars(coef))

# # Save the distribution plot
# ggsave(file.path(output_path,
#        "P_value_interactions_plot.pdf"),
#        plot = p_value_dist, height = 10, width = 6)

# # Filter significant interactions
# limma_res_sig <- limma_res |>
#   as_tibble() |>
#   filter(adj.P.Val < 0.05)

# # Plot feature with biggest interaction effect
# feature <- limma_res_sig |>
#   filter(str_detect(coef, ":")) |>
#   group_by(coef) |>
#   slice_max(logFC) |>
#   pull(Feature)

# # Get statistic
# # annotation <- limma_res_sig |>
# #   filter(str_detect(coef, ":")) |>
# #   group_by(coef) |>
# #   slice_max(logFC) |>
# #   select(coef, Feature, logFC , adj.P.Val) |>
# #   mutate(coef = gsub(":Breast_Invasive_Lobular_Carcinoma", "", coef)) |>
# #   mutate({{ancestry_column}} := coef)

# # Plot this one feature
# feature_expression <- as_tibble(t(data_norm$E[feature,]), rownames = 'sampleId')
# meta_data <- as_tibble(data$obs, rownames = 'sampleId')
# feature_meta <- merge(meta_data, feature_expression, on = 'sampleId')
# # Pivot longer
# feature_meta <- feature_meta |>
#   pivot_longer(cols = all_of(feature),
#                names_to = "Feature",
#                values_to = "Value")

# # Plot
# single_feature_plot <- feature_meta |>
#   ggplot(
#     aes(
#       x = get(output_column),
#       y = Value
#     )
#   ) +
#   geom_boxplot() +
#   facet_grid(cols = vars(get(ancestry_column)),
#              rows = vars(Feature)) +
#   xlab("Comparsion") + 
#   ylab("Normalized Expression") +
#   theme(axis.text.x = element_text(angle = 90))

# ggsave(file.path(output_path,
#        "Single_features.pdf"),
#        plot = p_value_dist, height = 10, width = 6)


# # Visualize multiple features
# # Get up reg. features
# goi_up <- limma_res_sig |>
#   filter(str_detect(coef, ":"), 
#          adj.P.Val < 0.05) |>
#   group_by(coef) |>
#   filter(logFC > 0) |>
#   slice_max(logFC, n = 5) |>
#   pull(Feature)

# # Get down reg. features
# goi_down <- limma_res_sig |>
#   filter(str_detect(coef, ":"), 
#          adj.P.Val < 0.05) |>
#   group_by(coef) |>
#   filter(logFC < 0) |>
#   slice_min(logFC, n = 5) |>
#   pull(Feature)

# # Statistic results
# # Plot up reg. features
# stat_up <- limma_res_interactions |>
#   filter(Feature %in% goi_up) |>
#   ggplot(
#     aes(
#       y = Feature,
#       x = str_replace(coef, ":.*$", ""),
#       color = logFC,
#       size = -log10(adj.P.Val))
#   ) + 
#   geom_point() +
#   scale_color_gradient2(high="red", low="blue") +
#   xlab("Comparison") +
#   ylab("Feature") +
#   ggtitle("Top upreg. interactions") + 
#   theme(axis.text.x = element_text(angle = 90))

# # Expression
# exp_list <- list()
# for(gg in goi_up){
#   exp_list[[gg]] <- data$obs |>
#     mutate(E=scale(data_norm$E[gg,])) |>
#     rownames_to_column("sampleId") |>
#     remove_rownames()
# }
# # Bind to dataframe
# expression_up <- bind_rows(exp_list, .id = "Feature")

# # Plot expression heatmap
# heatmap_up <- expression_up |>
#   group_by(cancer_type_detailed, genetic_ancestry, Feature) |>
#   summarise(across(where(is.numeric), mean)) |>
#   ggplot(
#     aes(
#       x = genetic_ancestry,
#       y = Feature,
#       fill = E
#     )
#   ) +
#   geom_tile() +
#   facet_grid(cols = vars(cancer_type_detailed)) +
#   scale_fill_gradient2(low="blue", high="red") +
#   ggtitle("Average expression per ancestry") 

# # Patchwork (combine stats results and expression)
# upregulated_plot <- heatmap_up + stat_up
# # Save the patchwork plot
# ggsave(file.path(output_path, 
#        "Upregulated_interactions.pdf"),
#        plot = upregulated_plot, height = 5, width = 10)

# # Plot down reg. features
# stat_down <- limma_res_interactions |>
#   filter(Feature %in% goi_down) |>
#   ggplot(
#     aes(
#       y = Feature,
#       x = str_replace(coef, ":.*$", ""),
#       color = logFC,
#       size = -log10(adj.P.Val))
#   ) + 
#   geom_point() +
#   scale_color_gradient2(high="red", low="blue") +
#   xlab("Comparison") +
#   ylab("Feature") +
#   ggtitle("Top downreg. interactions") + 
#   theme(axis.text.x = element_text(angle = 90))

# # Expression
# exp_list <- list()
# for(gg in goi_down){
#   exp_list[[gg]] <- data$obs |>
#     mutate(E=scale(data_norm$E[gg,])) |>
#     rownames_to_column("sampleId") |>
#     remove_rownames()
# }
# # Bind to dataframe
# expression_down <- bind_rows(exp_list, .id = "Feature")

# # Plot expression heatmap
# heatmap_down <- expression_down |>
#   group_by(cancer_type_detailed, genetic_ancestry, Feature) |>
#   summarise(across(where(is.numeric), mean)) |>
#   ggplot(
#     aes(
#       x = genetic_ancestry,
#       y = Feature,
#       fill = E
#     )
#   ) +
#   geom_tile() +
#   facet_grid(cols = vars(cancer_type_detailed)) +
#   scale_fill_gradient2(low="blue", high="red") +
#   ggtitle("Average expression per ancestry") 

# # Patchwork (combine stats results and expression)
# downregulated_plot <- heatmap_down + stat_down
# # Save the patchwork plot
# ggsave(file.path(output_path, 
#        "Downregulated_interactions.pdf"),
#        plot = downregulated_plot, height = 5, width = 10)









