# Remove start up messages
suppressPackageStartupMessages(
  {
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # DGE workflow and functional analysis
    library(edgeR)
    library(Rtsne)
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
    library(tidyverse)
    library(data.table)
    library(yaml)
    library(anndata)
  }
)
# Custom functions
source("r_utils.R")
source("figure_themes.R")

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

  # data/inputs/settings/PanCanAtlas_BRCA_RSEM_Basal_vs_non-Basal_EUR_to_ADMIX.yml
  # data/inputs/settings/PanCanAtlas_BRCA_RPPA_Basal_vs_non-Basal_EUR_to_ADMIX.yml
  # data/inputs/settings/Firehose_BRCA_BETA_Basal_vs_non-Basal_EUR_to_ADMIX.yml

  # data/inputs/settings/PanCanAtlas_LUSC_LUAD_RSEM_LUSC_vs_LUAD_EUR_to_ADMIX.yml
  # data/inputs/settings/PanCanAtlas_LUSC_LUAD_RPPA_LUSC_vs_LUAD_EUR_to_ADMIX.yml
  # data/inputs/settings/Firehose_LUSC_LUAD_BETA_LUSC_vs_LUAD_EUR_to_ADMIX.yml


  yaml_file <- "data/inputs/settings/PanCanAtlas_LUSC_LUAD_RPPA_LUSC_vs_LUAD_EUR_to_ADMIX.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}
# Set seed 
set.seed(42)

# Construction:
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out = "data/combined_runs"
analysis_name = "descriptive_statisitics"


# Transform settings into R useable form
tag             <- setup$tag
comparison      <- setup$classification$comparison
output_column   <- setup$classification$output_column
ancestry_column <- setup$classification$ancestry_column
train_ancestry  <- setup$classification$train_ancestry
inf_ancestry    <- setup$classification$infer_ancestry

# Construct output directory name
condition <- paste(comparison, collapse = "_vs_")

# Create directory if it does not exists
match_pattern <- paste0(tag, "_", condition, "_", analysis_name)
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Make directory
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Load data
adata <- read_h5ad(setup$data_path)

# Clustering 
# Transform data
if (setup$data_type == "protein") {

  cluster_data <- adata$X

  # Visualize
  p_before <- plot_density_of_samples(adata$X, x_axis_label = "Input values") 
  p_after <- plot_density_of_samples(cluster_data, x_axis_label = "Input values")
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_name <- file.path(path_to_save_location, "QC_density_values.pdf")
  save_ggplot(p, save_name, width = 3, height = 3)

} else if (setup$data_type == "expression") {
  # Expression (raw RSEM): LogCPM
  norm_factors <- calculate_tmm_norm_factors(adata$X)
  cluster_data <- cpm(adata$X, norm_factors = norm_factors, log = TRUE)

  # Visualize
  p_before <- plot_density_of_samples(adata$X, x_axis_label = "RSEM values") 
  p_after <- plot_density_of_samples(cluster_data, x_axis_label = "logCPM")
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_name <- file.path(path_to_save_location, "QC_density_values.pdf")
  save_ggplot(p, save_name, width = 3, height = 3)

} else if (setup$data_type == "methylation") {
  # Methylation (beta-values): M-values
  cluster_data <- beta_to_mvalue(adata$X)

  # Visualize
  p_before <- plot_density_of_samples(adata$X, x_axis_label = "Beta values") 
  p_after <- plot_density_of_samples(cluster_data, x_axis_label = "M-values")
  # Combine
  p <- p_before + p_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # Save
  save_name <- file.path(path_to_save_location, "QC_density_values.pdf")
  save_ggplot(p, save_name, width = 3, height = 3)

} else{
  cluster_data <- adata$X
}


# TSNE
set.seed(42)
# Reduce 
if (ncol(cluster_data) < 300){
  perplexity = 10
} else{
  perplexity = 50
}

tsne <- Rtsne(cluster_data, dims = 2, perplexity)
# Coordinates
tsne_coordinates <- data.frame(TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2])
# Add meta data
tsne_coordinates  <- bind_cols(tsne_coordinates, adata$obs)

# Visualize 
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
  theme_nature_fonts() +
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
  theme_nature_fonts() +
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
  theme_nature_fonts() +
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
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
  ) +
  guides(color = guide_legend(ncol = 1))

# ---- Patchwork ----
tsne_plot <- tsne_genetic_ancestry_plot +
  tsne_cancer_type_detailed_plot +
  tsne_subtype_plot +
  tsne_sex_plot +
  plot_layout(
    nrow = 2, ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.spacing = unit(2, "mm"), 
    legend.key.height = unit(5, "point"),         
    legend.key.width = unit(5, "point"),
    legend.margin = margin(0, 0, 0, 0)                   
  )

# Save
ggsave(filename = "Patchwork_clustering_plot.pdf",
       path = path_to_save_location,
       plot = tsne_plot,
       height = 3, width = 4.5)

# Save without legend
cluster_title <- str_to_sentence(setup$data_type)
tsne_plot_no_legend <- tsne_plot + 
  plot_annotation(title = cluster_title) & 
  theme(
    plot.title = element_text(hjust = 0.5, size = 5),
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
  theme_nature_fonts() +
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
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = c(legend_x_pos, legend_y_pos),  
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
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(legend.position = "none") 

# Save
ggsave(filename = "Variance_across_ancestry.pdf",
       plot = variance_across_ancestry_plot,
       path = path_to_save_location,
       height = 3, width = 3)
