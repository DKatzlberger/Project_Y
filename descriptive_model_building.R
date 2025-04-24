# Remove start up messages
suppressPackageStartupMessages(
  { 
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # Clustering
    library(Rtsne)
    # Parallelization
    library(parallel)
    library(furrr)
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
  YAML_FILE <- "example_settings_descriptive_model_building.yaml"
  is_yaml_file(YAML_FILE)
}

# Input
# --- Settings
# Default
DEFAULT_FILE  <- "default_settings_descriptive_model_building.yaml"
default_setup <- load_default_settings(DEFAULT_FILE)
# User
user_setup <- yaml.load_file(YAML_FILE)
# Default and user settings (user will overwrite default)
setup      <- modifyList(default_setup, user_setup)
# Check required settings
# Required settings
required_settings <- c(
  "data_path", "tech", "output_directory",
  "meta_variables"
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
# Load data
data_path  <- setup$data_path
is_h5ad_file(data_path)
adata      <- read_h5ad(setup$data_path)
# Check if columns exist
required_columns <- setup$meta_variables
check_columns(adata$obs, required_columns)

# --- Descriptive statistics
# Message
cat(sprintf("New analysis with id: %s; created: %s \n", setup$id, setup$date))
cat(sprintf("Save location: %s \n", path_to_save_location))
cat("-------------------------------------------------------------------- \n")

# --- Normalization/Transformation
tech      <- setup$tech
data_type <- setup$data_type

save_name <- file.path(path_to_save_location, "QC_density_values.pdf")
if (tech == "transcriptomics"){

  # Print statement
  cat(sprintf("Transforamtion for technology: %s. \n", tech))

  # LogCPM transformation
  norm_factors <- calculate_tmm_norm_factors(adata$X)
  trans_data   <- cpm(adata$X, norm_factors = norm_factors, log = TRUE)
  # Plot
  p_before <- plot_density_of_samples(adata$X, x_axis_label = data_type)
  p_after  <- plot_density_of_samples(trans_data, x_axis_label = "logCPM")
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

  # Mean variance trend
  # Axis
  x_axis <- paste0("log2(", data_type, " + 0.5)")
  # Data
  data_ <- log2(adata$X + 0.5)
  # Plot
  p <- plot_mean_variance_trend(data_, x_axis)
  # Save
  save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (tech == "methylation"){

  # Beta to mvals
  trans_data <- beta_to_mvalue(adata$X)

  # Plot
  p_before <- plot_density_of_samples(adata$X, x_axis_label = data_type)
  p_after  <- plot_density_of_samples(trans_data, x_axis_label = "M-values")
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

  # Mean variance trend
  # Axis
  x_axis <- paste0("log2(", data_type, " + 0.5)")
  # Data
  data_ <- log2(adata$X + 0.5)
  # Plot
  p <- plot_mean_variance_trend(data_, x_axis)
  # Save
  save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

} else if (tech == "proteomics"){

  # No transformation
  trans_data <- adata$X

  # Plot
  p_before <- plot_density_of_samples(adata$X, x_axis_label = data_type)
  p_after  <- plot_density_of_samples(trans_data, x_axis_label = data_type)
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_ggplot(p, save_name, width = 6, height = 4)

  # Mean variance trend
  # Axis
  x_axis <- paste0("log2(", data_type, " + 0.5)")
  # Data
  data_ <- log2(adata$X + 0.5)
  # Plot
  p <- plot_mean_variance_trend(data_, x_axis)
  # Save
  save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
  save_ggplot(p, save_name, width = 6, height = 4)

}
# Number of features
setup$n_features  <- ncol(trans_data)
cat("-------------------------------------------------------------------- \n")

# --- Unsupervised clustering 
cat("Clustering analysis. \n")
# Meta variables to analyse
meta_variables <- setup$meta_variables
# TSNE
cat("Unsuperived clustering: TSNE. \n")
set.seed(setup$seed)
tsne_results <- Rtsne(trans_data, dims = 2, setup$perplexity)
coordinates  <- data.frame(TSNE1 = tsne_results$Y[, 1], TSNE2 = tsne_results$Y[, 2])
coordinates  <- bind_cols(coordinates, adata$obs)
# Save
if (setup$save_coordinates){
  save_name <- file.path(path_to_save_location, "TSNE_coordinates.csv")
  fwrite(coordinates, save_name)
}
# Visualize: TSNE
p_list <- lapply(meta_variables, function(var) plot_clusters(coordinates, "TSNE1", "TSNE2", var, var))
# Patchwork
p <- wrap_plots(p_list)
# Save
save_name <- file.path(path_to_save_location, "TSNE_cluster.pdf")
save_ggplot(p, save_name, width = 5, height = 5)

# UMAP
cat("Unsuperived clustering: UMAP. \n")
library(umap) # Placed here because of error
umap_results <- umap(trans_data, n_components = 2, n_neighbors = setup$perplexity, seed = setup$seed)
coordinates  <- data.frame(UMAP1 = umap_results$layout[, 1],  UMAP2 = umap_results$layout[, 2])
coordinates  <- bind_cols(coordinates, adata$obs)
# Save 
if (setup$save_coordinates){
  save_name <- file.path(path_to_save_location, "UMAP_coordinates.csv")
  fwrite(coordinates, save_name)
}
# Visualize: UMAP
p_list <- lapply(meta_variables, function(var) plot_clusters(coordinates, "UMAP1", "UMAP2", var, var))
# Patchwork
p <- wrap_plots(p_list)
# Save
save_name <- file.path(path_to_save_location, "UMAP_cluster.pdf")
save_ggplot(p, save_name, width = 5, height = 5)

# PCA
cat("Unsuperived clustering: PCA. \n")
pca_results <- prcomp(trans_data, center = TRUE)
coordinates <- data.frame(PCA1 = pca_results$x[, 1],  PCA2 = pca_results$x[, 2])
coordinates <- bind_cols(coordinates, adata$obs)
# Save 
if (setup$save_coordinates){
  save_name <- file.path(path_to_save_location, "PCA_coordinates.csv")
  fwrite(coordinates, save_name)
}
# Visualize: PCA
p_list <- lapply(meta_variables, function(var) plot_clusters(coordinates, "PCA1", "PCA2", var, var))
# Patchwork
p <- wrap_plots(p_list)
# Save the combined PCA plot
save_name <- file.path(path_to_save_location, "PCA_cluster.pdf")
save_ggplot(p, save_name, width = 5, height = 5)
cat("-------------------------------------------------------------------- \n")

# --- Variance across PCs
cat("Variance explained by PCs. \n")
# Settings
meta_variables <- setup$meta_variables
n_pcs          <- setup$n_pcs
# Remove variable with less than 2 unique values
var_count    <- check_unique_values(adata$obs, meta_variables)
filtered_var <- var_count[var_count$Count >= 2, ]$Variable
# Combine PCs with meta 
meta_ <- as_tibble(adata$obs, rownames = "Idx")
meta_ <- select(meta_, Idx, all_of(filtered_var))
data_ <- as_tibble(pca_results$x[, 1:n_pcs], rownames = "Idx")
# Join and pivot 
pc_data <- inner_join(meta_, data_, by = "Idx")
pc_data <- pivot_longer(
  data      = pc_data,
  cols      = -c(Idx, all_of(filtered_var)),
  names_to  = "PC",
  values_to = "Expression"
)

# Lm (explain variance)
n_variables <- length(filtered_var)
n_features  <- n_pcs
batch_size  <- 10000
# Calculate number of cores
request_cores <- n_variables * ceiling(n_features / batch_size)
# Validate availabe cores
total_cores <- parallel::detectCores()
if (total_cores <= request_cores){
  # Save one core
  available_cores <- parallel::detectCores() - 1
} else {
  available_cores <- request_cores
}
# Function call
pcs_explained <- LM_explain(
  data       = pc_data, 
  variables  = filtered_var, 
  feature    = "PC",
  response   = "Expression",
  batch_size = batch_size,
  n_cores    = available_cores
)
# Add numeric column
pcs_explained <- mutate(
  pcs_explained, 
  feature_numeric = as.numeric(gsub("PC", "", feature)),
  feature_numeric = as.factor(feature_numeric)
)
# Save 
save_name <- file.path(path_to_save_location, "PCA_explained.csv")
fwrite(pcs_explained, save_name)

# Visualize: Variance across PCs
# Pivot
pcs_explained <- pivot_longer(
  data      = pcs_explained,
  cols      = c(R_squared, p_value),
  names_to  = "metric",
  values_to = "values"
)
# Loop for plots
p_list <- lapply(filtered_var, function(var) {
  # Filter the data
  data_r2 <- filter(pcs_explained, variable == var, metric == "R_squared")
  data_p  <- filter(pcs_explained, variable == var, metric == "p_value")
  # Plot
  plot_r2_pval(data_r2, data_p, x = "feature_numeric", title = var)
})
# Patchwork
p <- wrap_plots(p_list)
# Save the combined PCA plot
save_name <- file.path(path_to_save_location, "PCA_explained.pdf")
save_ggplot(p, save_name, width = 15, height = 10)


# --- Variance across features
# Merge meta with expression
meta_ <- as_tibble(adata$obs, rownames = "Idx")
meta_ <- select(meta_, Idx, all_of(filtered_var))
data_ <- as_tibble(trans_data, rownames = "Idx")
# Join and pivot
exp_data <- inner_join(meta_, data_, by = "Idx")
exp_data <- pivot_longer(
  data      = exp_data,
  cols      = -c(Idx, all_of(filtered_var)),
  names_to  = "Feature",
  values_to = "Expression"
)
# Lm (explain variance)
n_variables <- length(filtered_var)
n_features  <- length(unique(exp_data$Feature))
batch_size  <- 10000
# Calculate number of cores
request_cores <- n_variables * ceiling(n_features / batch_size)
# Validate availabe cores
total_cores <- parallel::detectCores()
if (total_cores <= request_cores){
  # Save one core
  available_cores <- parallel::detectCores() - 1
} else {
  available_cores <- request_cores
}
# Function call
features_explained <- LM_explain(
  data       = exp_data, 
  variables  = filtered_var, 
  feature    = "Feature",
  response   = "Expression",
  batch_size = batch_size,
  n_cores    = available_cores
)



# Save the settings
save_name <- file.path(path_to_save_location, "Settings.yaml")
write_yaml(setup, save_name)



print("Unsupervised clustering.")
point_size <-  0.5

# ---- Genetic ancestry (tsne_genetic_ancestry_plot) ----
genetic_ancestry_colors <- c(
  "admix" = "#ff4d4d",
  "afr"   = "#ff9900", 
  "amr"   =  "#33cc33",
  "eur"   = "#3399ff", 
  "eas"   = "#cc33ff", 
  "sas"   = "#ffcc00"
  )

tsne_genetic_ancestry_plot <- tsne_coordinates |>
  ggplot(
    aes(
      x = TSNE1,
      y = TSNE2,
      color = genetic_ancestry,
    )
  ) +
  geom_point(size = point_size) +
  scale_color_manual(
    values = genetic_ancestry_colors,
    labels = function(x) toupper(x)
  ) +
  labs(
    title = "Genetic ancestry",
    color = "Genetic ancestry"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
  ) +
  guides(color = guide_legend(ncol = 2))

# ---- Cancer type detailed (tsne_cancer_type_detailed_plot) ----
cancer_type_detailed_colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", 
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
  "#f1c40f", "#e74c3c", "#9b59b6", "#2ecc71"
)

tsne_cancer_type_detailed_plot <- tsne_coordinates |>
  ggplot(
    aes(
      x = TSNE1,
      y = TSNE2,
      color = cancer_type_detailed
    )
  ) +
  geom_point(size = point_size) +
  scale_color_manual(
    values = cancer_type_detailed_colors
  ) +
  labs(
    title = "Histological subtype",
    color = "Histological subtype"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
  ) +
  guides(color = guide_legend(ncol = 1))

# ---- Subtype (tsne_subtype_plot) ----
subtype_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", 
  "#66a61e", "#e6ab02", "#a6761d", "#666666"
)
tsne_subtype_plot <- tsne_coordinates |>
  ggplot(
    aes(
      x = TSNE1,
      y = TSNE2,
      color = subtype
    )
  ) +
  geom_point(size = point_size) +
  scale_color_manual(
    values = subtype_colors
  ) +
  labs(
    title = "Molecular subtype",
    color = "Molecular subtpye"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(ncol = 2))

# ---- Sex (tsne_sex_plot) ----
tsne_sex_plot <- tsne_coordinates |>
  ggplot(
    aes(
      x = TSNE1,
      y = TSNE2,
      color = sex
    )
  ) +
  geom_point(size = point_size) +
  scale_color_manual(
    values = c("Female" = "#c7640d", "Male" = "#0e67a7")
  ) +
  labs(
    title = "Sex",
    color = "Sex"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
  ) +
  guides(color = guide_legend(ncol = 1))

# ---- Patchwork ----
cluster_title <- str_to_sentence(setup$data_type)
# Plot
tsne_plot <- tsne_genetic_ancestry_plot +
  tsne_cancer_type_detailed_plot +
  tsne_subtype_plot +
  tsne_sex_plot +
  plot_annotation(
    title = cluster_title
  ) +
  plot_layout(
    nrow = 2, 
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    legend.spacing = unit(2, "mm"), 
    legend.key.height = unit(5, "point"),         
    legend.key.width = unit(5, "point"),
    legend.margin = margin(0, 0, 0, 0)                   
  )

# Save
ggsave(filename = "Patchwork_clustering_plot.pdf",
       path = path_to_save_location,
       plot = tsne_plot,
       height = 3.5, width = 4.5)

# Save without legend
tsne_plot_no_legend <- tsne_plot + 
  plot_annotation(
    title = cluster_title
  ) & 
  theme(
    plot.title = element_text(hjust = 0.5, size = font_size),
    legend.position = "none"
    )

ggsave(filename = "Patchwork_clustering_plot_no_legend.pdf",
       path = path_to_save_location,
       plot = tsne_plot_no_legend,
       height = 3, width = 3)


# Define classification task
adata <- adata[adata$obs[[output_column]] %in% comparison]

# Descriptive statistics 
print("Descriptive analsyis.")
alpha_value <- 0.3

# ---- Genetic ancestry (genetic_ancestry_plot) -----
genetic_ancestry_plot <- adata$obs |>
  ggplot(
    aes(
      x = !!sym(output_column),
      fill = !!sym(ancestry_column)
    )
  ) +
  geom_bar(position = "fill") +
  coord_flip() + 
  scale_y_continuous(
    breaks = c(0.0, 0.5, 1.0),
    labels = c(0, 0.5, 1)
  ) +
  scale_fill_manual(
    values = genetic_ancestry_colors,
    labels = function(x) toupper(x)
  ) +
  facet_grid(
    rows = vars(!!sym(ancestry_column))
  ) +
  labs(
    title = "Genetic ancestry",
    x = "Condition",
    y = "Proportion",
    fill = "Genetic ancestry"
  ) +
  theme_nature_fonts() + 
  theme_small_legend() + 
  theme_white_background() + 
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    strip.text = element_blank(),      
    strip.background = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  ) +
  guides(fill = guide_legend(ncol = 3))

# ---- Cancer type detailed (cancer_type_detailed_plot) ----
cancer_type_detailed_colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", 
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
)
cancer_type_detailed_plot <- adata$obs |>
  ggplot(
    aes(
      x = !!sym(setup$classification$output_column),
      fill = cancer_type_detailed
    )
  ) +
  geom_bar(position = "fill") +
  coord_flip() + 
  scale_y_continuous(
    breaks = c(0.0, 0.5, 1.0),
    labels = c(0, 0.5, 1)
  ) +
  scale_fill_manual(
    values = cancer_type_detailed_colors
  ) +
  facet_grid(
    rows = vars(!!sym(setup$classification$ancestry_column))
  ) +
  labs(
    title = "Histological\nsubtype",
    x = "Condition",
    y = "Proportion",
    fill = "Histological subtype"
  ) +
  theme_nature_fonts() + 
  theme_small_legend() + 
  theme_white_background() + 
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    strip.text = element_blank(),      
    strip.background = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  ) +
  guides(fill = guide_legend(ncol = 1))

# ---- Subtype (subtype_plot) ----
subtype_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", 
  "#66a61e", "#e6ab02", "#a6761d", "#666666"
)
subtype_plot <- adata$obs |> 
  ggplot(
    aes(
      x = !!sym(setup$classification$output_column),
      fill = subtype
    )
  ) +
  geom_bar(position = "fill") +
  coord_flip() + 
  scale_y_continuous(
    breaks = c(0.0, 0.5, 1.0),
    labels = c(0, 0.5, 1) 
  ) +
  scale_fill_manual(values = subtype_colors) +
  facet_grid(
    rows = vars(!!sym(setup$classification$ancestry_column))
    ) +
  labs(
    title = "Molecluar\nsubtype",
    x = "Comparison",
    y = "Proportion",
    fill = "Molecluar subtype"
  ) +
  theme_nature_fonts(
    # base_size = 5
  ) + 
  theme_small_legend() + 
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    strip.text = element_blank(),      
    strip.background = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  ) +
  guides(fill = guide_legend(ncol = 2))

# ---- Sex (sex_plot) -----
sex_plot <- adata$obs |> 
  ggplot(
    aes(
      x = !!sym(setup$classification$output_column),
      fill = sex
    )
  ) +
  geom_bar(
    position = "fill"
    ) +
  coord_flip() + 
  scale_y_continuous(
    breaks = c(0.0, 0.5, 1.0),
    labels = c(0, 0.5, 1)
  ) +
  scale_fill_manual(
    values = c("Female" = "#c7640d", "Male" = "#0e67a7")
  ) +
  facet_grid(
    rows = vars(!!sym(setup$classification$ancestry_column))
    ) +
  labs(
    title = "Sex",
    x = "Comparison",
    y = "Proportion",
    fill = "Sex"
  ) +
  theme_nature_fonts(
    #base_size = 5
  ) +
  theme_small_legend() + 
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    strip.text = element_blank(),      
    strip.background = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  ) +
  guides(fill = guide_legend(ncol = 3))

# ---- Total plot (total_plot) -----
total_plot <- adata$obs |> 
  ggplot(
    aes(
      x = !!sym(setup$classification$output_column),
    )
  ) +
  geom_bar() +
  geom_text(
    stat = "count", 
    aes(label = after_stat(count)), 
    hjust = -0.3, 
    size = 1.5
  ) +
  coord_flip() + 
  scale_y_continuous(
    n.breaks = 3,
    expand = expansion(mult = c(0.1, 0.35))
  ) +
  facet_grid(rows = vars(!!sym(setup$classification$ancestry_column))) +
  labs(
    title = "Total patients",
    x = "Comparison",
    y = "Total",
  ) +
  theme_nature_fonts() + 
  theme_white_background() +
  theme(
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    strip.text = element_blank(),      
    strip.background = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  )

# ---- Age plot (age_plot) ----
age_plot <- adata$obs |> 
  ggplot(
    aes(
      x = age,
      fill = !!sym(setup$classification$output_column)
    )
  ) +
  geom_density(
    alpha = 0.5,
    color = NA
    ) +
  scale_y_continuous(
    n.breaks = 3
  ) +
  facet_grid(
    rows = vars(toupper(!!sym(setup$classification$ancestry_column))),
    scale = "free"
  ) +
  labs(
    title = "Age", 
    x = "Age",
    y = "Density",
    fill = "Condition"
  ) +
  theme_nature_fonts() +
  theme_small_legend() + 
  theme_white_background() +
  theme_white_strip() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    axis.text.y = element_blank(),  
    axis.title.y = element_blank(),
    plot.margin = margin(5, 1, 5, 1),
  ) +
  guides(fill = guide_legend(ncol = 1))

# ---- Patchwork ----
nf_plot <- genetic_ancestry_plot +
  cancer_type_detailed_plot + 
  subtype_plot + 
  sex_plot + 
  total_plot + 
  age_plot + 
  plot_layout(ncol = 6, guides = "collect") &
  theme(
    legend.spacing = unit(2, "mm"), 
    legend.key.height = unit(5, "point"),         
    legend.key.width = unit(5, "point"),
    legend.margin = margin(0, 0, 0, 0)                      
  ) 

# Save
ggsave(filename = "Patchwork_descriptive_meta_plot.pdf",
       path = path_to_save_location,
       plot = nf_plot,
       height = 3, width = 9)

# Save without the legend
nf_plot_no_legend <- nf_plot & theme(legend.position = "none")

# Save the plot without the legend
ggsave(filename = "Patchwork_descriptive_meta_plot_no_legend.pdf",
       path = path_to_save_location,
       plot = nf_plot_no_legend,
       height = 3, width = 3)


# Variance 
print("Variance explained.")
# Variance by meta variable
check_variance <- c(
  ancestry_column,
  "subtype",
  "cancer_type_detailed", 
  "subtype",
  "sex", 
  "age"
  )

# Meta data
filtered_meta <- adata$obs |>
  as_tibble(rownames = "Idx") |>
  dplyr::select(Idx, all_of(check_variance))

# Expression
filtered_expression <- adata$X |>
  as_tibble(rownames = "Idx") |>
  inner_join(filtered_meta, by = "Idx")

# Reshape
filtered_expression <- filtered_expression |>
  pivot_longer(
    -c(Idx, all_of(check_variance)),
    names_to = "Feature",
    values_to = "Expression"
  )

# R2 
variance_explained <- calculate_variance_explained(
  expression_df = filtered_expression,
  check_variance = check_variance,
  batch_size = 1000,
  num_cores = setup$njobs
) 

# Additional information
variance_explained <- variance_explained |>
  mutate(
    data_type = setup$data_type
  )

# Save
fwrite(variance_explained, file.path(path_to_save_location, "Variance_explained.csv"))


# Variance by ancestry
variance_by_ancestry <- data.table()
for (condition in unique(adata$obs[[output_column]])){
  # Meta data
  filtered_meta <- adata$obs |>
    as_tibble(rownames = "Idx") |>
    filter(!!sym(output_column) == condition) |>
    dplyr::select(Idx, all_of(ancestry_column))

  # Expression
  filtered_expression <- adata$X |>
    as_tibble(rownames = "Idx") |>
    inner_join(filtered_meta, by = "Idx")

  # Reshape
  filtered_expression <- filtered_expression |>
    pivot_longer(
      -c(Idx, all_of(ancestry_column)),
      names_to = "Feature",
      values_to = "Expression"
    )

  # R2 per condition
  con_variance_by_ancestry <- calculate_variance_explained(
    expression_df = filtered_expression,
    check_variance = c(ancestry_column),
    batch_size = 1000,
    num_cores = setup$njobs
  ) |>
  mutate(Condition = condition)

  # Combine
  variance_by_ancestry <- bind_rows(variance_by_ancestry, con_variance_by_ancestry)
}

# Additional information
variance_by_ancestry <- variance_by_ancestry |>
  mutate(
    data_type = setup$data_type
  )

# Save
fwrite(variance_by_ancestry, file.path(path_to_save_location, "Variance_by_ancestry.csv"))

# Variance across ancestries
variance_across_ancestry <- data.table()
for (condition in unique(adata$obs[[output_column]])){
  # Meta data
  filtered_meta <- adata$obs |>
    as_tibble(rownames = "Idx") |>
    filter(!!sym(output_column) == condition) |>
    dplyr::select(Idx, all_of(ancestry_column))
  
  # Expression
  filtered_expression <- adata$X |>
    as_tibble(rownames = "Idx") |>
    inner_join(filtered_meta, by = "Idx")
  
  # Reshape
  filtered_expression <- filtered_expression |>
    pivot_longer(
      -c(Idx, all_of(ancestry_column)),
      names_to = "Feature",
      values_to = "Expression"
    )

  # Calculate variance
  con_variance_across_ancestry <- filtered_expression |>
    group_by(!!sym(ancestry_column), Feature) |>
    summarize(
      Variance = var(Expression, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(Condition = condition)
  
  # Combine
  variance_across_ancestry <- bind_rows(variance_across_ancestry, con_variance_across_ancestry)
}


# Visualize 
# Mapping
variance_mapping <- c(
  "genetic_ancestry"      = "Genetic ancestry",
  "cancer_type_detailed"  = "Histologocal subtype",
  "subtype"               = "Molecular subtype",
  "sex"                   = "Sex",
  "age"                   = "Age"
)

variance_colors <- c(
  "Genetic ancestry"      = "#1f77b4", 
  "Histologocal subtype"  = "#ff7f0e",  
  "Molecular subtype"     = "#2ca02c",  
  "Sex"                   = "#d62728",  
  "Age"                   = "#9467bd"  
)

variance_explained <- variance_explained %>%
  mutate(
    mapped_variable = recode(Variable, !!!variance_mapping),
  )

# Legend x-position
max_name_length <- max(nchar(unique(variance_explained$mapped_variable)))
legend_x_pos <- max(0.5, 0.9 - 0.005 * max_name_length)
# Legend y-position
max_items <- length(unique(variance_explained$mapped_variable))
legend_y_pos <- max(0.3, 0.9 - 0.01 * max_items)

# Relevel variables
variance_explained <- variance_explained %>%
  mutate(
    mapped_variable = factor(
      mapped_variable, 
      levels = c(
        "Genetic ancestry", 
        "Histologocal subtype", 
        "Molecular subtype", 
        "Sex", 
        "Age"
      )
    ) 
  )

# Summarize
avg_r2_var <- variance_explained[, .(Mean_R2 = mean(R2, na.rm = TRUE)), by = .(Variable)]
avg_r2_cond <- variance_by_ancestry[, .(Mean_R2 = mean(R2, na.rm = TRUE)), by = .(Variable = Condition)]
# Merge the two datasets
avg_r2 <- rbind(avg_r2_var, avg_r2_cond, use.names = TRUE, fill = TRUE)


# ---- Variance explained (variance_explained_plot) ----
variance_explained_plot <- variance_explained |>
  filter(!is.na(R2)) |>
  ggplot(
    aes(
      x = R2,
      fill = mapped_variable
    )
  ) +
  geom_histogram(
    binwidth = 0.001, 
    alpha = alpha_value,
    position = "identity"
  ) +
  scale_x_continuous(
    limits = c(0.0, 0.5),  
    breaks = seq(0, 0.5, by = 0.1)  
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = variance_colors
  ) +
  labs(
    x = "Proportion of variance explained",
    y = "Number of genes",
    fill = "Variable"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = c(legend_x_pos, legend_y_pos),  
    legend.direction = "vertical",
    legend.title = element_blank()
  ) 

# Save
ggsave(filename = "Variance_explained.pdf",
       plot = variance_explained_plot,
       path = path_to_save_location,
       height = 1.5, width = 3)

# ---- Variance by ancestry (variance_by_ancestry_plot) ----
variance_by_ancestry_plot <- variance_by_ancestry |>
  left_join(
    avg_r2_cond,
    by = c("Condition" = "Variable")
  ) |>
  mutate(
    Condition = paste0(Condition, " ~ ", round(Mean_R2 * 100, 2))
  ) |>
  ggplot(
    aes(
      x = R2,
      fill = Condition
    )
  ) +
  geom_histogram(
    binwidth = 0.001, 
    alpha = 0.7,
    position = "identity"
  ) +
  scale_x_continuous(
    limits = c(0.0, 0.5),  
    breaks = seq(0, 0.5, by = 0.1)  
  ) +
  scale_fill_manual(
    values = c("lightblue", "lightpink")
  ) +
  labs(
    x = "Proportion of variance explained",
    y = "Number of genes",
    fill = "Condition"
  ) +
  theme_nature_fonts(font_size) +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = c(0.8, 0.8),  
    legend.direction = "vertical",
    legend.title = element_blank()
  ) 

# Save
ggsave(filename = "Variance_by_ancestry.pdf",
       plot = variance_by_ancestry_plot,
       path = path_to_save_location,
       height = 1.5, width = 3)

# ---- Variance across ancestry (variance_across_ancestry_plot) ----
variance_across_ancestry_plot <- variance_across_ancestry |>
  filter(Variance > 0) |>
  ggplot(
    aes(
      x = toupper(!!sym(ancestry_column)),
      y = Variance
    )
  ) +
  geom_violin(
    aes(fill = !!sym(ancestry_column)), 
    alpha = alpha_value
  ) +
  geom_boxplot(
    width = 0.2, 
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = genetic_ancestry_colors
  ) +
  scale_y_log10() +
  facet_grid(
    rows = vars(Condition),
    scales = "free"
  ) +
  labs(
    x = "Ancestry",
    y = "Variance"
  ) +
  theme_nature_fonts(font_size) +
  theme_white_background() +
  theme_white_strip() +
  theme(legend.position = "none") 

# Save
ggsave(filename = "Variance_across_ancestry.pdf",
       plot = variance_across_ancestry_plot,
       path = path_to_save_location,
       height = 3, width = 3)
