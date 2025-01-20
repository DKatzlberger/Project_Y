# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Statistics 
    library(coin)
    # Visualization
    library(patchwork)
    }
)
# Source custom functions
source("r_utils.R")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  yaml_file <- args[1] 
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
  # Vscratch directory
  vscratch_dir_in = file.path("data", "runs")
  # When the script is run from the command line then 'output_directory' is given
  # The pattern to extract all matchin directories is extracted from 'output_directory'
  output_path = setup$output_directory
  match_pattern <- sub("_\\d+$", "", sub(".*/", "", output_path))
} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "PanCanAtlas_RSEM_covariates_subtype_age_EUR_to_EAS.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

  # Construction:
  # Vscratch_dir is the place where the files are stored
  vscratch_dir_in = file.path("data", "runs")
  # Tag is used to specify which data it is e.g. TCGA, NIAGADS
  tag <- setup$tag
  # Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
  comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
  # Analysis_name specifies what was analysed
  # E.g. comparsion of ancestries EUR_to_AMR, subsetting_EUR, robustness_AFR
  # This often is modified depending which analysis
  train_ancestry <- toupper(setup$classification$train_ancestry)
  infer_ancestry <- toupper(setup$classification$infer_ancestry)
  analysis_name  <- paste0(train_ancestry, "_to_", infer_ancestry, "_robustness")
  # Combination of components to create the 'match_pattern'
  # The 'match_pattern' is used as pattern to extract all folders in the vscratch dir
  match_pattern <- paste0(comparison, "_", analysis_name)
}

# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)
# Print the match folders
print("Matched folders:")
print(match_vscratch_dir)

# Check if there were matching folders
if (length(match_vscratch_dir) == 0) {
  message("No matching folders found.")
}

# Save the results of the analysis
# 'vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# Combining metric csv files 
# Each folder is one seed
metric_dge <- data.frame()
metric_ml <- data.frame()
for (folder in match_vscratch_dir){
    dge_file <- file.path(folder, "Metric_dge.csv")
    ml_file <- file.path(folder, "Metric_ml.csv")

    # Load and append DGE data for each seed
    dge_data <- fread(dge_file) 
    metric_dge <- bind_rows(metric_dge, dge_data) 

    # Load and append ML data for each seed
    ml_data <- fread(ml_file) 
    metric_ml <- bind_rows(metric_ml, ml_data) 
}

# Summarize metric
# 1. Permutation testing
# 2. Correct for multiple testing
summarized_dge <- tibble()
summarized_ml <- tibble()
for (prop in setup$proportion) {
  # Filter by proportion
  filtered_dge <- filter(metric_dge, Proportion == prop)
  filtered_ml  <- filter(metric_ml, Proportion == prop)

  # Perform permutation testing (non-parametric)
  # dge
  p_pearson <- permutation_test(data = filtered_dge, group_col = "Status", value_col = "Pearson")
  p_spearman <- permutation_test(data = filtered_dge, group_col = "Status", value_col = "Spearman")
  # ml
  p_regression <- permutation_test(data = filtered_ml, group_col = "Status", value_col = "LogisticRegression")
  p_forest <- permutation_test(data = filtered_ml, group_col = "Status", value_col = "RandomForestClassifier")

  # Summarize 'filtered_dge'
  summarized_filtered_dge <- filtered_dge |>
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = "Value",
        names_to = "Metric"
    ) |>
    group_by(
      Ancestry, 
      Status,
      Prediction, 
      Metric, 
      n_inf_ancestry, 
      Proportion
    ) |>  
    summarize(
      n_seeds = n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      se_value = sd(Value, na.rm = TRUE) / sqrt(n())
    ) |>
    # Add p_value
    mutate(p_value = ifelse(Metric == "Pearson", p_pearson, p_spearman))

    # Summarize 'filtered_ml'
    summarized_filtered_ml <- filtered_ml |>
      pivot_longer(
        cols = c(LogisticRegression, RandomForestClassifier),
        values_to = "Value",
        names_to = "Algorithm"
    ) |> 
      group_by(
        Ancestry, 
        Status,
        Prediction, 
        Algorithm, 
        n_inf_ancestry, 
        Proportion
      ) |>  
      summarize(
        n_seeds = n(),
        mean_value = mean(Value, na.rm = TRUE),
        sd_value = sd(Value, na.rm = TRUE),
        se_value = sd(Value, na.rm = TRUE) / sqrt(n())
      ) |>
      # Add p_value
      mutate(p_value = ifelse(Algorithm == "LogisticRegression", p_regression, p_forest))

    # Combine data
    summarized_dge <- bind_rows(summarized_dge, summarized_filtered_dge)
    summarized_ml <- bind_rows(summarized_ml, summarized_filtered_ml)
}


# Pvalue correction:
# Correct for multiple testing (Benjamini-Hochberg correction)
# dge 
unique_p_values_dge <- summarized_dge |> pull(p_value) |> unique()
p_adjusted_bh_dge <- p.adjust(unique_p_values_dge, method = "BH")
# Create a named vector to map the adjusted p-values back to the data frame
p_adjusted_named_dge <- setNames(p_adjusted_bh_dge, unique_p_values_dge)
# Add to 'summarized_dge'
summarized_dge <- summarized_dge |>
   mutate(p_adjusted = p_adjusted_named_dge[as.character(p_value)])

# ml 
unique_p_values_ml <- summarized_ml |> pull(p_value) |> unique()
p_adjusted_bh_ml <- p.adjust(unique_p_values_ml, method = "BH")
# Create a named vector to map the adjusted p-values back to the data frame
p_adjusted_named_ml <- setNames(p_adjusted_bh_ml, unique_p_values_ml)
# Add to 'summarized_ml'
summarized_ml <- summarized_ml |>
   mutate(p_adjusted = p_adjusted_named_ml[as.character(p_value)])

# Save the metric data frames
fwrite(metric_dge, file.path(path_to_save_location, "Metric_dge.csv"))
fwrite(metric_ml, file.path(path_to_save_location, "Metric_ml.csv"))
# 2. Summarized metric
fwrite(summarized_dge, file.path(path_to_save_location, "Summarized_metric_dge.csv"))
fwrite(summarized_ml, file.path(path_to_save_location, "Summarized_metric_ml.csv"))


# Visulaization:
# dge
# 1. Create a plot for each proportion
# 2. Pathwork to one plot
ggplot_list_dge = list()
for (prop in setup$proportion) {
  # Filter by proportion
  filtered_dge <- summarized_dge |> filter(Proportion == prop)

  # Create lable to showcase number of samples
  n_inf_label <- filtered_dge |>
    ungroup() |>  
    distinct(Ancestry, n_inf_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
    select(Ancestry, label) |>
    deframe()
  
  # Correlation axis from 0 to 1
  common_y <- scale_y_continuous(
    limits = c(0, 1.1), 
    breaks = c(0, 0.5, 1))
  
  # Create plot
  prop_plot <- filtered_dge |>
    ggplot(
      aes(
        x = fct_rev(Prediction),
        y = mean_value
        )
    ) +
    geom_bar(
      stat = "identity", 
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        ymin = mean_value - sd_value, 
        ymax = mean_value + sd_value
        ),
      width = 0.2, 
      position = position_dodge(0.7)
    ) +
    common_y +
    facet_grid(
      rows = vars(Metric),
      col = vars(Ancestry),
      labeller = labeller(Ancestry = as_labeller(n_inf_label), 
                          Metric = label_value)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "Correlation coefficient"
    ) +
    geom_text(
      aes(
        x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
        y = Inf,  # Position the text at the top of the plot
        label = paste("Perm. test,", "p =", format(p_adjusted, digits = 3))
        ),
      size = 4,    
      vjust = 2,   # Align text to the top (since we're using Inf for y)
      hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
      inherit.aes = FALSE  
  ) 

  # Add the plot to the list with a unique name
  ggplot_list_dge[[paste0("prop_", prop)]] <- prop_plot
}

# Combine with patchwork 
combined_dge_plot <- wrap_plots(ggplot_list_dge, ncol = 2)
# Save the image
ggsave(filename = "Plot_dge.pdf", 
       plot = combined_dge_plot, 
       path = path_to_save_location, 
       width = 10, height = 10
       )

# ml 
ggplot_list_ml = list()
for (prop in setup$proportion) {
  # Filter by proportion
  filtered_ml <- summarized_ml |> filter(Proportion == prop)

  # Create lable to showcase number of samples
  n_inf_label <- filtered_ml |>
    ungroup() |>  
    distinct(Ancestry, n_inf_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ")")) |>
    select(Ancestry, label) |>
    deframe()
  
  # Correlation axis from 0 to 1
  common_y <- scale_y_continuous(
    limits = c(0, 1.1), 
    breaks = c(0, 0.5, 1))
  
  # Create plot
  prop_plot <- filtered_ml |>
    ggplot(
      aes(
        x = fct_rev(Prediction),
        y = mean_value
        )
    ) +
    geom_bar(
      stat = "identity", 
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        ymin = mean_value - sd_value, 
        ymax = mean_value + sd_value
        ),
      width = 0.2, 
      position = position_dodge(0.7)
    ) +
    common_y +
    facet_grid(
      rows = vars(Algorithm),
      col = vars(Ancestry),
      labeller = labeller(Ancestry = as_labeller(n_inf_label), 
                          Metric = label_value)
    ) +
    labs(
        # title = paste(gsub("_", " ", comparison)),
        x = "Prediction",
        y = "ROC AUC"
    ) +
    geom_text(
      aes(
        x = 0.5,  # Align text to the left side of the first bar (0.5 is a bit to the left of the first bar)
        y = Inf,  # Position the text at the top of the plot
        label = paste("Perm. test,", "p =", format(p_adjusted, digits = 3))
        ),
      size = 4,    
      vjust = 2,   # Align text to the top (since we're using Inf for y)
      hjust = 0,   # Align text to the left (since we're using x = 0.5 for the first bar)
      inherit.aes = FALSE  
  ) 

  # Add the plot to the list with a unique name
  ggplot_list_ml[[paste0("prop_", prop)]] <- prop_plot
}

# Combine with patchwork 
combined_ml_plot <- wrap_plots(ggplot_list_ml, ncol = 2)
# Save the image
ggsave(filename = "Plot_ml.pdf", 
       plot = combined_ml_plot, 
       path = path_to_save_location, 
       width = 10, height = 10
       )


