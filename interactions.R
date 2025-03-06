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
    library(Rtsne)
    library(clusterProfiler)
    library(org.Hs.eg.db)
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
  yaml_file <- "data/inputs/settings/PanCanAtlas_BRCA_RPPA_basal_vs_non-basal_EUR_to_AFR.yml"
  print("Running interactive mode for development.")
  setup <- yaml.load_file(yaml_file)
}
# Set seed 
set.seed(42)

# Construction:
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out = "data/combined_runs"
analysis_name = "interactions"

tag <- setup$tag
comparison <- paste(setup$classification$comparison, collapse = "_vs_")
train_ancestry <- toupper(setup$classification$train_ancestry)
inf_ancestry <- toupper(setup$classification$infer_ancestry)
ancestry <- paste0(train_ancestry, "_to_", inf_ancestry)

# Match pattern
match_pattern <- paste0(tag, "_", comparison, "_", ancestry, "_", analysis_name)
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Transform settings into R useable form
comparison <- setup$classification$comparison
output_column <- setup$classification$output_column
ancestry_column <- setup$classification$ancestry_column
train_ancestry <- setup$classification$train_ancestry
inf_ancestry <- setup$classification$infer_ancestry

# Load the data
adata_whole <- read_h5ad(setup$data_path)

# Clustering -------------------------------------------------------------------
# Transform data to log2 
if (setup$data_type == "rppa") {
  # For RPPA (raw data), apply shift to handle negative values
  min_value <- min(adata_whole$X, na.rm = TRUE)
  shifted_rppa <- adata_whole$X + abs(min_value) + 1 
  cluster_data <- log2(shifted_rppa)  
} else if (setup$data_type == "expression") {
  # Log2 normalization
  cluster_data <-  log2(cluster_data + 1)
} else {
  # Default: Keep input data as is
  cluster_data <- adata_whole$X
}
# Clustering (tsne)
tsne <- Rtsne(cluster_data, dims = 2, perplexity = 50)
# Coordinates
tsne_coordinates <- data.frame(TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2])
# Add meta data
tsne_coordinates  <- bind_cols(tsne_coordinates, adata_whole$obs)

# Visualize 
point_size <-  0.5

# ---- Genetic ancestry (tsne_genetic_ancestry_plot) ----
genetic_ancestry_colors <- c("admix"= "#ff4d4d", "afr" = "#ff9900", "amr" =  "#33cc33",
                             "eur" = "#3399ff", "eas" = "#cc33ff", "sas" = "#ffcc00"
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
    values = genetic_ancestry_colors
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
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
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
    title = "Molecular-subtype",
    color = "Molecular-subtpye"
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
tsne_plot_no_legend <- tsne_plot & theme(legend.position = "none")

ggsave(filename = "Patchwork_clustering_plot_no_legend.pdf",
       path = path_to_save_location,
       plot = tsne_plot_no_legend,
       height = 3, width = 3)


# Define classification task
adata <- adata_whole[adata_whole$obs[[output_column]] %in% comparison]


# Descriptive statistics 
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

# --- Patchwork ----
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
# Variance by meta variable
check_variance <- c(
  ancestry_column, 
  "cancer_type_detailed", "subtype",
  "sex", "age"
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
  "genetic_ancestry" = "Genetic ancestry",
  "cancer_type_detailed" = "Histologocal subtype",
  "subtype" = "Molecular subtype",
  "sex" = "Sex",
  "age" = "Age"
)

variance_colors <- c(
  "Genetic ancestry" = "#1f77b4", 
  "Histologocal subtype" = "#ff7f0e",  
  "Molecular subtype" = "#2ca02c",  
  "Sex" = "#d62728",  
  "Age" = "#9467bd"  
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
    mapped_variable = factor(mapped_variable, levels = c("Genetic ancestry", 
                                                        "Histologocal subtype", 
                                                        "Molecular subtype", 
                                                        "Sex", 
                                                        "Age")) 
  )


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
    limits = c(0.0, 0.3),  
    breaks = seq(0, 0.3, by = 0.1)  
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
       height = 3, width = 3)

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
       height = 3, width = 3)

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


# Differential gene expression --------------------------------------------------
# The intercept in the 'Interaction' analysis should be the same as in 'cross_ancestry'
# Setting intercept of ancestry: EUR
# Setting intercept of comparison: Based on occurance in 'setup$classification$comparison'

# Filtering by ancestry
data <- adata[adata$obs[[ancestry_column]] %in% c(train_ancestry, inf_ancestry)]

# Relevel 'ancestry_column' to be 'train_ancestry' as the intercept
# 1. Get possible values for ancestries
# 2. Relevel them so 'train_ancestry' is first level
# 3. Update in the actual data frame
possible_ancestry <- unique(data$obs[[ancestry_column]])
releveled_ancestry <- factor(
  possible_ancestry, 
  levels = c(train_ancestry, setdiff(possible_ancestry, train_ancestry))
  )
levels_ancestry <- levels(releveled_ancestry)

# Update 'ancestry_column' in the actual data
data$obs <- data$obs |>
  mutate(
    !!ancestry_column := factor(
      .data[[ancestry_column]],
      levels = levels_ancestry)
      )

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

# Store baseline 
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

# Rename to R useable names
colnames(design) <- make.names(colnames(design))

# Filter features 
if (setup$filter_features & setup$data_type == "expression"){
  # Filter by expression
  keeper_genes <- filterByExpr(t(data$X), design = design)
  filtered_data <- data[, keeper_genes]
} else{
  filtered_data <- data
}
# Save feature 
write_yaml(filtered_data$var_names, file.path(path_to_save_location, "Features.yml"))

# Limma workflow
print("Start differential gene expression analysis.")
# Select normalization method
data_type <- setup$data_type
dge_normalization <- setup$dge_normalization
normalization_method <- normalization_methods[[data_type]][[dge_normalization]]
# Transpose (rows = Genes, cols = Samples)
filtered_data <- t(filtered_data$X)
# Normalization
norm_data <- normalization_method(filtered_data, train_design)

# Fit the model (means model)
limma_fit <- lmFit(norm_data, design = design)
limma_fit <- eBayes(limma_fit)
limma_fit_res <- extract_results(limma_fit)
# Remove interecpt
limma_fit_res <- filter(limma_fit_res, coef != "X.Intercept.")

# Baseline --------------------------------------------------------------------------------------------------
# Condition which is used for the intercept
intercept_baseline <- limma_fit_res |> 
  filter(coef == inf_ancestry) |>
  mutate(coef = str_replace_all(coef, inf_ancestry, "intercept_baseline")) |>
  mutate(Condition = setup$classification$comparison[1])

# Non-intercept baseline (need to fit contrast)
contrast_term <- paste(inf_ancestry, "+", paste0(inf_ancestry, ".", setup$classification$comparison[2]))
contrast.matrix <- makeContrasts(
  contrast_term,
  levels = design
)
colnames(contrast.matrix) <- "non_intercept_baseline"

# Fit contrast
limma_fit_contrast <- contrasts.fit(limma_fit, contrast.matrix)
limma_fit_contrast <- eBayes(limma_fit_contrast)
# Extract results
non_intercept_baseline <- extract_results(limma_fit_contrast) |>
  mutate(Condition = setup$classification$comparison[2])

# Combine baseline comparisons
baseline <- bind_rows(intercept_baseline, non_intercept_baseline) |>
  mutate(
    Ancestry = toupper(inf_ancestry),
    Comparison = paste(toupper(train_ancestry), "to", toupper(inf_ancestry))
  )
# Save 
fwrite(baseline, file.path(path_to_save_location, "Baseline.csv"))

# Define significant baseline differences
logFC_threshold <- 1
sig_baseline <- baseline |>
  filter(
    adj.P.Val < 0.05 & 
    abs(logFC) > logFC_threshold
  ) |>
  pull(Feature)

# Interactions --------------------------------------------------------------------------------------------------
interactions <- filter(limma_fit_res, str_detect(coef, "\\.")) |>
  mutate(
    Ancestry = toupper(inf_ancestry),
    Comparison = paste(toupper(train_ancestry), "to", toupper(inf_ancestry)),
    Condition = "Interactions"
  )
# Save the interactions
fwrite(interactions, file.path(path_to_save_location, "Interactions.csv"))

# Visualize 
logFC_threshold <- 1
alpha_value  <- 0.7
point_size  <- 0.5

# Sig. interactions
sig_interactions <- interactions |>
  filter(
    adj.P.Val < 0.05 & 
    abs(logFC) > logFC_threshold
  ) |>
  pull(Feature)

# Top interactions
top_interactions <- interactions |>
  filter(
    adj.P.Val < 0.05 & 
    abs(logFC) > logFC_threshold
  ) |>
  slice_max(abs(logFC), n = 5) |>
  pull(Feature)

# ---- Baseline volcano plot (baseline_volcano_plot) ----
baseline_volcano_plot <- baseline |>
  ggplot(
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
    )
  ) +
  geom_point(
    size = point_size
  ) +
  geom_text_repel(
    data = baseline |> filter(Feature %in% top_interactions),
    aes(label = Feature),
    size = 1.5,
    color = "black",  
    segment.color = "black",  
    min.segment.length = 0
  ) +
  facet_grid(
    cols = vars(Ancestry),
    rows = vars(Condition)
  ) + 
  geom_vline(
    xintercept = c(-logFC_threshold, logFC_threshold), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = alpha_value
  ) +
  geom_hline(
    yintercept = -log10(0.05), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = alpha_value
  ) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "lightgrey"),
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Baseline_volcano_plot.pdf",
       plot = baseline_volcano_plot,
       path = path_to_save_location,
       height = 3, width = 3)

# ---- Interactions volcano plot (interactions_volcano_plot) ----
interactions_volcano_plot <- interactions |>
  ggplot(
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = (adj.P.Val < 0.05 & abs(logFC) > logFC_threshold)
    )
  ) +
  geom_point(
    size = point_size
  ) +
  geom_text_repel(
    data = interactions |> filter(Feature %in% top_interactions),
    aes(label = Feature),
    size = 1.5,
    color = "black",  
    segment.color = "black",  
    min.segment.length = 0
  ) +
  facet_grid(
    cols = vars(Ancestry),
    rows = vars(Condition)
  ) + 
  geom_vline(
    xintercept = c(-logFC_threshold, logFC_threshold), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = alpha_value
  ) +
  geom_hline(
    yintercept = -log10(0.05), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = alpha_value
  ) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "lightgrey"),
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_white_strip() +
  theme(legend.position = "none")

# Save 
ggsave(filename = "Interactions_volcano_plot.pdf",
       plot = interactions_volcano_plot,
       path = path_to_save_location,
       height = 3, width = 3)

# ---- Correlation of baselines (baseline_scatter_plot) -------
# Reformat data
wide_baseline <- baseline |> 
  mutate(Comparison = paste(Comparison, Condition)) |>
  dplyr::select(Feature, Comparison, logFC) |>
  pivot_wider(names_from = Comparison, values_from = logFC) |>
  rename_with(~ str_replace_all(., " ", "_"))

# Fit linear model and get predictions
x_col <- names(wide_baseline)[2]
y_col <- names(wide_baseline)[3]
lm_model <- lm(as.formula(paste(y_col, "~", x_col)), data = wide_baseline)

# Coefficients
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]

predicted_values <- fitted(lm_model)

residuals <- residuals(lm_model)
abs_residuals <- abs(residuals)

# Define outliers
residual_threshold <- 1
upper_threshold <- predicted_values + residual_threshold
lower_threshold <- predicted_values - residual_threshold

# Calculate the correlation
cor_pearson <- cor(wide_baseline[[x_col]], wide_baseline[[y_col]], method = "pearson")
cor_spearman <- cor(wide_baseline[[x_col]], wide_baseline[[y_col]], method = "spearman")
# Label
label_pearson <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
label_spearman <- bquote(R[Spearman] == .(round(cor_pearson, 3)))

# Reformat data
scatter_plot_data <- baseline |>
  mutate(Comparison = paste(Comparison, Condition, "(logFC)")) |>
  dplyr::select(Feature, Comparison, logFC) |>
  pivot_wider(names_from = Comparison, values_from = logFC) |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > residual_threshold,
    upper_threshold = upper_threshold,
    lower_threshold = lower_threshold
  ) |>
  mutate(
    with_interactions = Feature %in% sig_interactions, 
    with_baseline = Feature %in% sig_baseline,
    coloring = case_when(
      with_interactions & with_baseline ~ "Interaction + Baseline",
      with_interactions &! with_baseline ~ "Interaction",    
      with_baseline & !with_interactions ~ "Baseline",
      TRUE ~ "Rest"                                
    )
  )

# Residuals above threshold
sig_residuals <- scatter_plot_data |>
  filter(abs(residuals) > logFC_threshold) |>
  pull(Feature)

# Legend 
# Legend x-position
max_name_length <- max(nchar(unique(scatter_plot_data$coloring)))
legend_x_pos <- min(0.5, 0.05 + 0.005 * max_name_length)
# Legend y-position
max_items <- length(unique(scatter_plot_data$coloring))
legend_y_pos <- max(0.3, 0.9 - 0.01 * max_items)

# Plot
baseline_scatter_plot <- scatter_plot_data %>%
  ggplot(
    aes(
      x = .data[[names(.)[2]]],
      y = .data[[names(.)[3]]],
      color = coloring
    )
  ) +
  geom_point(
    data = scatter_plot_data |> filter(coloring == "Rest"),
    aes(color = coloring),
    size = point_size,
    show.legend = FALSE
  ) +
  geom_point(
    data = scatter_plot_data |> filter(coloring == "Rest"),
    aes(color = "dummy"),
    size = NA,   
    show.legend = FALSE 
  ) +
  geom_point(
    data = scatter_plot_data |> filter(coloring == "Baseline"),
    aes(color = coloring),
    size = point_size
  ) +
  geom_point(
    data = scatter_plot_data |> filter(coloring %in% c("Interaction + Baseline", "Interaction")),
    aes(color = coloring),
    size = point_size
  ) +
  geom_smooth(
    method = "lm", 
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) + 
  geom_line(
    aes(y = upper_threshold), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +  
  geom_line(
    aes(y = lower_threshold), 
    linetype = "dashed", 
    color = "blue",
    linewidth = point_size,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Interaction + Baseline" = "red",
      "Interaction" = "darkgreen", 
      "Baseline" = "gold", 
      "Rest" = "lightgrey",
      "dummy" = "black"
    ),
    breaks = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      "Rest",
      "dummy"
    ),
    labels = c(
      "Interaction + Baseline", 
      "Interaction", 
      "Baseline",
      label_pearson,
      label_spearman
    )
  ) +
  scale_x_continuous(
    limits = c(
      min(scatter_plot_data[[names(scatter_plot_data)[2]]]) - 1, 
      max(scatter_plot_data[[names(scatter_plot_data)[2]]]) + 1
    )  
  ) +
  scale_y_continuous(
    limits = c(
      min(scatter_plot_data[[names(scatter_plot_data)[3]]]) - 1, 
      max(scatter_plot_data[[names(scatter_plot_data)[3]]]) + 1
    )  
  ) +
  geom_text_repel(
    data = scatter_plot_data |> filter(Feature %in% top_interactions),
    aes(label = Feature),
    size = 1.5,
    color = "black",  
    segment.color = "black",  
    min.segment.length = 0
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend() +
  theme(
    legend.position = c(legend_x_pos, legend_y_pos), 
    legend.title = element_blank(),  
    legend.spacing.x = unit(0, "cm"), 
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_rect(fill = "transparent"),  
    legend.key.spacing = unit(0, "cm")
  )

# Save the plot
ggsave(filename = "Baseline_scatter_plot.pdf",
       plot = baseline_scatter_plot,
       path = path_to_save_location,
       height = 3, width = 3)

# ---- Venn Diagram (venn_diagram) ----
venn_data <- list(
  Baseline = sig_baseline,
  Interactions = sig_interactions,
  Residuals = sig_residuals
)

venn_data <- Venn(venn_data)
venn_data <- process_data(venn_data)

# Adjust coordinates of labels
venn_data_adjusted <- venn_setlabel(venn_data) %>%
  mutate(
    adjusted_X = case_when(
      X == max(X) ~ X * 0.8,  #
      X == min(X) ~ X * 0.8,  #
      TRUE ~ X  
    ),
    adjusted_Y = case_when(
      X == max(X) ~ Y * 0.8,  
      X == min(X) ~ Y * 1.2,  
      TRUE ~ Y 
    )
  )

# Now plot with adjusted X and Y values
vennDiagram <- ggplot() +
  # 1. region count layer
  geom_polygon(
    data = venn_regionedge(venn_data),
    aes(X, Y, group = id),
    fill = "white"
  ) +
  # 2. set edge layer
  geom_path(
    data = venn_setedge(venn_data), 
    aes(X, Y, group = id), 
    show.legend = FALSE
  ) +
  # 3. set label layer
  geom_text(
    data = venn_data_adjusted,
    aes(adjusted_X, Y, label = name),
    size = 2
  ) +
  # 4. region label layer
  geom_label(
    data = venn_regionlabel(venn_data),
    aes(X , Y, label = count), 
    size = 2
  ) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

# Save
ggsave(filename = "Venn_diagram.pdf",
       plot = vennDiagram,
       path = path_to_save_location,
       height = 3, width = 3)

# ---- Heatmap (expression_heatmap) ----
exp_list <- list()
for(gg in sig_interactions){
  exp_list[[gg]] <- data$obs |>
    mutate(E=scale(data_norm$E[gg,])) |>
    rownames_to_column("patient_id") |>
    remove_rownames()
}
# Bind to dataframe
expression <- bind_rows(exp_list, .id = "Feature")

pdf(file.path(path_to_save_location, "Interactions_expression_heatmap.pdf"), 
    width = 6 , height = 3
  )
# Expression matrix
expression_matrix <- expression |>
  dplyr::select(Feature, patient_id, E) |>
  spread(key = patient_id, value = E) |>
  column_to_rownames("Feature") |>
  as.matrix()

# Annotations
# Colors
output_column_colors <- setNames(c("#1b9e77", "purple"), comparison)

# Annotation dataframe
anno_df <- expression |>
  dplyr::select(patient_id, all_of(c(output_column, ancestry_column))) |>
  distinct() |>
  arrange(match(patient_id, colnames(expression_matrix)))

# Grouping 
custom_order <- c(
  paste(train_ancestry, comparison[1], sep = "_"),  
  paste(train_ancestry, comparison[2], sep = "_"),  
  paste(setdiff(anno_df[[ancestry_column]], train_ancestry), comparison[1], sep = "_"),  
  paste(setdiff(anno_df[[ancestry_column]], train_ancestry), comparison[2], sep = "_")  
)

anno_df <- anno_df |>
  mutate(
    group = paste(anno_df[[ancestry_column]], anno_df[[output_column]], sep = "_"),
    group = factor(group, levels = custom_order)
  ) 

# Annotation (Sample)
top_anno <- HeatmapAnnotation(
  # Annotations
  Ancestry = anno_empty(
    border = FALSE, 
    height = unit(0.2, "cm"), 
    show_name = TRUE
  ),
  Condition = anno_df[[output_column]], 
  # Colors
  col = list(
    Condition = output_column_colors
  ), 
  # Fonts
  annotation_name_gp = gpar(fontsize = 5),
  border = TRUE,
  show_annotation_name = TRUE,
  # Height
  simple_anno_size = unit(0.2, "cm"),
  # Legend
  show_legend = FALSE
)

# Legends
ancetsry_lgd <- Legend(
  title = "Ancestry", 
  at = c(
    toupper(train_ancestry), 
    toupper(inf_ancestry)
    ), 
  legend_gp = gpar(
    fill = c(
      genetic_ancestry_colors[train_ancestry], 
      genetic_ancestry_colors[inf_ancestry]
      )
    ),
  # Font
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize = 5, fontface = "plain"),
  # Size
  grid_height = unit(0.3, "cm"), 
  grid_width = unit(0.3, "cm"),
  # Direction
  direction = "horizontal"
  )

condition_lgd <- Legend(
  title = "Condition", 
  at = c(
    comparison[1], 
    comparison[2]
    ), 
  legend_gp = gpar(
    fill = c(
      output_column_colors[comparison[1]], 
      output_column_colors[comparison[2]]
      )
    ),
  # Font
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize = 5, fontface = "plain"),
  # Size
  grid_height = unit(0.3, "cm"), 
  grid_width = unit(0.3, "cm"),
  # Direction
  direction = "horizontal"
  )

lgd_list <- c(ancetsry_lgd, condition_lgd)

# Heatmap
heatmap <- Heatmap(
  expression_matrix,
  name = "Expression", 
  show_column_names = FALSE,
  # Clustering
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  # Split
  column_split = anno_df$group,
  cluster_column_slices = FALSE,
  column_title = NULL,
  # Annotations
  top_annotation = top_anno,
  # Fonts
  row_names_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 5),
  # Legends
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5, fontface = "plain"), 
    labels_gp = gpar(fontsize = 5, fontface = "plain"),
    grid_height = unit(0.3, "cm"), 
    grid_width = unit(0.3, "cm"),
    direction = "horizontal"
    )
  )

# Add legend
draw(
  heatmap, 
  # Legend
  annotation_legend_list = lgd_list,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE
  )

# Spanning groups
group_block_anno = function(group, empty_anno, gp = gpar(), 
    label = NULL, label_gp = gpar()) {

    seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
    loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
    loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

    seekViewport("global")
    grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
        just = c("left", "bottom"), gp = gp)
    if(!is.null(label)) {
        grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
    }
}

# Add the blocks
group_block_anno(
  1:2, "Ancestry", 
  gp = gpar(fill = genetic_ancestry_colors[train_ancestry])
  )
group_block_anno(
  3:4, "Ancestry", 
  gp = gpar(fill = genetic_ancestry_colors[inf_ancestry])
  )
# Save off
dev.off()  

# ---- Boxplot interactions (interaction_boxplot) ----
exp_list <- list()
for(gg in top_interactions){
  exp_list[[gg]] <- data$obs |>
    mutate(E=scale(data_norm$E[gg,])) |>
    rownames_to_column("patient_id") |>
    remove_rownames()
}
# Bind to dataframe
expression <- bind_rows(exp_list, .id = "Feature")

interactions_boxplot_list <- list()
for (feature in unique(expression$Feature)){
  # Filter by feature
  filtered_expression <- expression |>
  filter(Feature == feature)

  feat_interactions_boxplot <- filtered_expression |>
    ggplot(
      aes(
        x = !!sym(output_column),
        y = E
      )
    ) +
    geom_boxplot(
      width = point_size,
      outlier.size = 0.1,
    ) +
    scale_y_continuous(
      breaks = scales::breaks_extended(n = 3)
    ) +
    facet_grid(
      cols = vars(fct_rev(toupper(!!sym(ancestry_column)))),
    ) +
    labs(
      title = feature,
      x = "Condition"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Append to list 
  interactions_boxplot_list[[feature]] <- feat_interactions_boxplot
}

# Patchwork
interactions_boxplot <- wrap_plots(
  interactions_boxplot_list,
  nrow = 1, 
  ncol = ceiling(length(interactions_boxplot_list) / 1)
)

# Save
ggsave(filename = "Interactions_boxplot.pdf",
       plot = interactions_boxplot,
       path = path_to_save_location,
       height = 3, width = 6)


# Functional analysis -----------------------------------------------------------------------------------------
# Fgsea
database <- "data/downloads/geneset-libraries/MSigDB_Hallmark_2020.txt"
database <- read_enrichR_database(database)
enrichment <- perform_gsea(interactions, database, rank_by = "logFC") |>
  mutate(Ancestry = toupper(inf_ancestry))
# Save
fwrite(enrichment, file.path(path_to_save_location, "Interactions_enrichment.csv"))

# Cluster profiler 
GO_annotations <- fread("data/downloads/geneset-libraries/GO_annotation.csv")
GO_ids <- inner_join(interactions, GO_annotations, by=c("Feature"="gene_name"))
# Filter significant genes
GO_sig_interactions <- filter(GO_ids, adj.P.Val < 0.05) |> pull(gene_id)
# Background (EUR condition1 vs condition2)
GO_background <- GO_ids |> pull(gene_id)
GO_background <- setdiff(GO_background, GO_sig_interactions)
# Remove input
GO_background <- setdiff(GO_background, GO_sig_interactions)

# Overrepresentation
GO_overrepresentation <- enrichGO(
  gene = GO_sig_interactions,
  universe = GO_background,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db, 
  ont = "ALL",
  pAdjustMethod = "BH", 
  qvalueCutoff = 0.05, 
  readable = TRUE)
GO_results  <- as_tibble(GO_overrepresentation) |>
  mutate(Ancestry = toupper(inf_ancestry))


# Visulaize ----------------------------------------------------------------------------------------------
# ---- Fgsea (interactions_fgsea_plot) -----
interactions_fgsea_plot <- enrichment |>
  arrange(desc(abs(NES))) |>
  filter(row_number() <= 15) |>
  ggplot(
    aes(
      x = Ancestry,
      y = pathway,
      color = NES,
      size = pmin(-log10(padj), 5)
    )
  ) +
  geom_point() +
  scale_size_binned(
     range = c(1, 3)    
  ) +
  scale_color_gradient2(
      high = "red", 
      mid = "white", 
      low = "blue"
  ) +
  labs(
    title = "Enrichment of interactions",
    y = "MSigDB Hallmark 2020 gene set",
    size = "-log10(adj.P.Val)"
  ) +
  theme_nature_fonts() +
  theme_small_legend() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title.position = "top",
    # legend.key.height = unit(0.4, "cm"),  
    # legend.key.width = unit(0.4, "cm"),
  )

# Save
ggsave(filename = "Interactions_fgsea_plot.pdf",
       plot = interactions_fgsea_plot,
       path = path_to_save_location,
       height = 3, width = 4)

# ---- GO terms (interactions_GO_term_plot) -----
interactions_GO_term_plot <- GO_results |>
  ggplot(
    aes(
      x = Ancestry,
      y = Description,
      size = -log(p.adjust),
      color = FoldEnrichment
      )
  ) +
  geom_point() +
  scale_size_binned(
     range = c(1, 3)    
  ) +
  scale_color_gradient2(
      high = "red", 
      mid = "white", 
      low = "blue"
  ) +
  labs(
    title = "GO overrepresentation",
    y = "GO terms",
    size = "-log10(adj.P.Val)"
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.key.height = unit(0.4, "cm"),  
    legend.key.width = unit(0.4, "cm"),
  )

# Save
ggsave(filename = "Interactions_GO_plot.pdf",
       plot = interactions_GO_term_plot,
       path = path_to_save_location,
       height = 3, width = 4)

# ---- Patchwork ----
interactions_patchwork <- interactions_volcano_plot +
  interactions_fgsea_plot +
  interactions_GO_term_plot +
  plot_layout(ncol = 2)

# Save
ggsave(filename = "Patchwork_interactions_plot.pdf",
       plot = interactions_patchwork,
       path = path_to_save_location,
       height = 6, width = 7)




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









