# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Parallelization
    library(parallel)
    # Statistics 
    library(coin)
    library(fgsea)
    # Visualization
    library(patchwork)
    library(ggrepel)
    }
)
# Source custom functions
source("r_utils.R")
source("figure_themes.R")


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
  yaml_file <- "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_EAS.yml"
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

# Save location
vscratch_dir_out  <- file.path("data", "combined_runs")
path_to_save_location <- file.path(vscratch_dir_out, match_pattern)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}

# ---------------------------------------------------------------------------------------------------------------------
# Extracting all folders in the 'vscratch_dir_in' that match 'match_pattern'
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vscratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Contrast
# DGE
contrast_metric_dge <- fload_data(match_vscratch_dir, file = "Contrast_metric_dge.csv")
raw_metric_dge <- fload_data(match_vscratch_dir, file = "LogFCs.csv")

# ML
contrast_metric_ml <- fload_data(match_vscratch_dir, file = "Contrast_metric_ml.csv")
raw_metric_ml <- fload_data(match_vscratch_dir, file = "Probabilities.csv")

# Statistics (for each proportion)
# DGE
p_pearson_list <- list()
p_spearman_list <- list()
for (proportion in unique(contrast_metric_dge$Proportion)){
    # Filter based on proportion
    filtered_contrast_metric_dge <- filter(contrast_metric_dge, Proportion == proportion)

    p_pearson <- permutation_test(
        filtered_contrast_metric_dge, 
        value = "Pearson", 
        group = "Status",
        paired = "Seed",
        tail = "left"
        )

    p_spearman <- permutation_test(
        filtered_contrast_metric_dge, 
        value = "Spearman", 
        group = "Status",
        paired = "Seed",
        tail = "left"
        )
    
    # Append the results to the lists with the proportion as the name
    p_pearson_list[[as.character(proportion)]] <- p_pearson
    p_spearman_list[[as.character(proportion)]] <- p_spearman
}
# Create dataframes
p_pearson_df <- data.frame(
  Proportion = as.numeric(names(p_pearson_list)),
  p_value = unlist(p_pearson_list),
  Metric = "Pearson"
)

p_spearman_df <- data.frame(
  Proportion = as.numeric(names(p_spearman_list)),
  p_value = unlist(p_spearman_list),
  Metric = "Spearman"
)

# Combine data frames
p_values_dge <- bind_rows(p_pearson_df, p_spearman_df)

# ML
p_regression_auc_list <- list()
p_regression_acc_list <- list()
p_forest_auc_list <- list()
p_forest_acc_list <- list()
for (proportion in unique(contrast_metric_ml$Proportion)){
    # Filter based on proportion
    filtered_contrast_metric_ml <- filter(contrast_metric_ml, Proportion == proportion)

    # Filter by algorithm
    regression <- filter(filtered_contrast_metric_ml, Algorithm == "LogisticRegression")
    forest <- filter(filtered_contrast_metric_ml, Algorithm == "RandomForestClassifier")

    # Logistic regression
    p_regression_auc <- regression |>
        permutation_test(
            value = "ROC_AUC", 
            group = "Status",
            paired = "Seed",
            tail = "left"
            )   

    p_regression_acc <- regression |>
        permutation_test(
            value = "Accuracy", 
            group = "Status",
            paired = "Seed",
            tail = "left"
            )
    
    # Random forest
    p_forest_auc <- forest |>
            permutation_test(
            value = "ROC_AUC", 
            group = "Status",
            paired = "Seed",
            tail = "left"
            )
    
    p_forest_acc <- forest |>
        permutation_test(
            value = "Accuracy", 
            group = "Status",
            paired = "Seed",
            tail = "left"
            )
    
    # Append the results to the lists with the name of the proportion
    p_regression_auc_list[[as.character(proportion)]] <- p_regression_auc
    p_regression_acc_list[[as.character(proportion)]] <- p_regression_acc
    p_forest_auc_list[[as.character(proportion)]] <- p_forest_auc
    p_forest_acc_list[[as.character(proportion)]] <- p_forest_acc
}
# Create dataframes
p_regression_auc_df <- data.frame(
  Proportion = as.numeric(names(p_regression_auc_list)),
  p_value = unlist(p_regression_auc_list),
  Metric = "ROC_AUC",
  Algorithm = "LogisticRegression"
)

p_regression_acc_df <- data.frame(
  Proportion = as.numeric(names(p_regression_acc_list)),
  p_value = unlist(p_regression_acc_list),
  Metric = "Accuracy",
  Algorithm = "LogisticRegression" 
)

p_forest_auc_df <- data.frame(
  Proportion = as.numeric(names(p_forest_auc_list)),
  p_value = unlist(p_forest_auc_list),
  Metric = "ROC_AUC",
  Algorithm = "RandomForestClassifier"
)

p_forest_acc_df <- data.frame(
  Proportion = as.numeric(names(p_forest_acc_list)),
  p_value = unlist(p_forest_acc_list),
  Metric = "Accuracy",
  Algorithm = "RandomForestClassifier"
)

# Combine data frames
p_values_ml <- bind_rows(p_regression_auc_df,
  p_regression_acc_df,
  p_forest_auc_df,
  p_forest_acc_df
  )

# Aggregate data
# DGE - contrast_metric_dge
aggregated_contrast_metric_dge <- contrast_metric_dge |>
  pivot_longer(
    cols = c(
      Pearson, 
      Spearman,
      # RMSE,
      # R2
      ),
    values_to = "Value",
    names_to = "Metric"
  ) |> 
  group_by(
    Proportion,
    Metric, 
    Ancestry, 
    Status, 
    Prediction, 
    n_inf_ancestry
  ) |>  
  summarize(
    n_seeds = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  left_join(p_values_dge, by = c("Metric", "Proportion"))

# DGE - raw_metric_dge
aggregated_raw_metric_dge <- raw_metric_dge |>
  filter(
    coef == unique(raw_metric_dge$coef)[1]
  ) |>
  group_by(
    Proportion,
    coef,
    Feature,
    Status,
    Ancestry,
    n_train_ancestry,
    n_inf_ancestry,
  ) |>
  summarise(
    n_seeds = n(),
    mean_value = mean(logFC, na.rm = TRUE),
    sd_value = sd(logFC, na.rm = TRUE),
    se_value = sd(logFC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ML - contrast_metric_ml
aggregated_contrast_metric_ml <- contrast_metric_ml |>
 pivot_longer(
    cols = c(
      ROC_AUC, 
      Accuracy
      ),
    values_to = "Value",
    names_to = "Metric"
  ) |> 
  group_by(
    Proportion,
    Metric,
    Algorithm,
    Ancestry, 
    Status, 
    Prediction, 
    n_inf_ancestry
  ) |>  
  summarize(
    n_seeds = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    se_value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  left_join(
    p_values_ml, by = c("Proportion", "Algorithm", "Metric")
  )

# Visualize
point_size <- 0.5
# Define limits for y-axis
max_y_value <- 1.0 # not considering error bars
max_sd <- 0.1 # estimated max sd
max_multiplicator <- 1.3
max_hight <- max_multiplicator + max_sd
# y-limit
max_y_limit <- max_y_value * max_multiplicator
# y-scale
common_y <- scale_y_continuous(
  limits = c(0.0, max_y_limit),  
  breaks = c(0.0, 0.5, 1.0)  
)
# Plots
# ---- Contrast correlation logFC (DGE) ----
gg_proportion_list <- list()
for (proportion in rev(unique(aggregated_contrast_metric_dge$Proportion))){
  # Filter by proportion
  filtered_data <- filter(aggregated_contrast_metric_dge, Proportion == proportion)
  # Create labels
  inf_n_label <- filtered_data |>
    distinct(Ancestry, Proportion, n_inf_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ", ", Proportion, ")")) |>
    select(Ancestry, label) |>
    deframe()

  # Y-value for p_value bar
  y_max <- max(
      filtered_data$mean_value + filtered_data$sd_value, 
      na.rm = TRUE
      )

  # Plot
  proportion_plot <- filtered_data |>
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
    geom_text(
      aes(
          label = round(mean_value, 3),  
          y = mean_value * 0.5  
      ),
      size = 2,  
      color = "white",  
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
      labeller = labeller(
        Ancestry = as_labeller(inf_n_label)
        )
    ) +
    labs(
      x = "Prediction",
      y = "Y"
    ) +
    annotate(
      "segment",
      x = 1, xend = 2,  
      y = y_max * 1.1, yend = y_max * 1.1,  
      linewidth = 0.3
    ) +
    annotate(
        "segment",
        x = 1, xend = 1,
        y = y_max * 1.08, yend = y_max * 1.105,
        linewidth = 0.3
    ) +
    annotate(
        "segment",
        x = 2, xend = 2,
        y = y_max * 1.08, yend = y_max * 1.105,
        linewidth = 0.3
    ) +
    geom_text(
        aes(
            x = 1.5,
            y = y_max * max_multiplicator,
            label = ifelse(
            is.na(p_value),
            "p[paired~perm.~test] == NA",
            paste0(
                "p[paired~perm.~test] == ", 
                ifelse(p_value < 0.01, 
                    format(p_value, digits = 3, scientific = TRUE), 
                    format(p_value, digits = 3)),
                ifelse(p_value < 0.01, " ~`**`",  # Two asterisks for p < 0.001
                    ifelse(p_value < 0.05, " ~`*`", ""))  # One asterisk for p < 0.05
            )
            )
        ),
        parse = TRUE,  
        size = 2
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip()

  # Append
  gg_proportion_list[[as.character(proportion)]] <- proportion_plot
}

# Patchwork
combined_proportions <- wrap_plots(gg_proportion_list, ncol = 2)

# Save
ggsave(filename = "Contrast_correlation_logFC.pdf", 
       plot = combined_proportions, 
       path = path_to_save_location, 
       width = 4, height = 5
       )

# ---- Contrast distribution avg. logFC (DGE) ----
gg_proportion_list <- list()
for (proportion in rev(unique(aggregated_raw_metric_dge$Proportion))){
  # Filter by proportion
  filtered_data <- filter(aggregated_raw_metric_dge, Proportion == proportion)

  # Train - Inference 
  # Create label
  inf_n_label <- filtered_data |>
    filter(Status == "Inference") |>
    distinct(Ancestry, Proportion, n_inf_ancestry) |>
    mutate(label = paste0(Ancestry, " (n = ", n_inf_ancestry, ", ", Proportion, ")")) |>
    select(Proportion, label) |>
    deframe()

  # Prepare data for plot
  lm_data <- filtered_data |>
    select(-c(Ancestry, sd_value, se_value)) |>
    pivot_wider(names_from = Status, values_from = mean_value)

  # Model
  lm_model <- lm(Inference ~ Train, data = lm_data)
  # Create threshold line
  thr = 1
  lm_data <- lm_data |>
    mutate(
      predicted = fitted(lm_model),  
      upper_threshold = predicted + thr,  
      lower_threshold = predicted - thr
    )

  # Correlation
  cor_pearson <- cor(lm_data$Train, lm_data$Inference, method = "pearson")
  cor_spearman <- cor(lm_data$Train, lm_data$Inference, method = "spearman")
  # P-value
  cor_p_pearson <- pvalue(independence_test(Inference ~ Train, data = lm_data))
  cor_p_spearman <- pvalue(independence_test(rank(Inference) ~ rank(Train), data = lm_data))
  # Dataframe
  inference <- data.frame(
    Pearson = cor_pearson, 
    Spearman = cor_spearman,
    Prediction = "Ancestry",
    Status = "Inference"
    )
  # Label
  pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
  spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))
  
  # Train and inf ancestry
  x_ancestry <- unique(filtered_data$Ancestry[filtered_data$Status == "Train"])
  y_ancestry <- unique(filtered_data$Ancestry[filtered_data$Status == "Inference"])

  # Compared conditions
  condition_1 <- unique(raw_metric_dge$coef)[1]
  condition_2 <- unique(raw_metric_dge$coef)[2]

  # Plot
  prop_train_inference <- lm_data |>
    ggplot(
      aes(
        x = Train,
        y = Inference
      )
    ) +
    geom_point(
      size = point_size,
      color = "lightgrey"
    ) +
    geom_point(
      aes(color = "dummy1"),
      size = NA
    ) +
    geom_point(
      aes(color = "dummy2"),
      size = NA
    ) +
    geom_smooth(
      method = "lm",
      color = "blue",
      linewidth = point_size,
      alpha = 0.5
    ) +
    geom_line(
      aes(x = Train, y = upper_threshold), 
      color = "blue", 
      linetype = "dashed",
      linewidth = point_size,
      alpha = 0.5
    ) + 
    geom_line(
      aes(x = Train, y = lower_threshold), 
      color = "blue", 
      linetype = "dashed",
      linewidth = point_size,
      alpha = 0.5
    ) + 
    scale_color_manual(
      values = c(
        "dummy1" = "black",
        "dummy2" = "black"
        ),
      breaks = c(
        "dummy1",
        "dummy2"
        ),
      labels = c(
        pearson_label,
        spearman_label
      )
    ) +
    guides(
      color = guide_legend(
        title = "Ancestry",
        override.aes = list(size = 0, alpha = 0)
      )
    ) +
    labs(
      x = paste(x_ancestry, condition_1, "vs", condition_2, "(avg. logFC)"),
      y = paste(y_ancestry, condition_1, "vs", condition_2, "(avg. logFC)")
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      legend.title = element_text(face = "bold"),
      legend.justification = c(0, 1),
      legend.position = c(0.05, 1.0), 
      legend.spacing.x = unit(0, "cm"), 
      legend.spacing.y = unit(0, "cm"),
      legend.background = element_rect(fill = "transparent"),  
      legend.key.spacing = unit(0, "cm")
    )


  # Train - Test

  # Prepare data for plot
  lm_data <- filtered_data |>
    select(-c(Ancestry, sd_value, se_value)) |>
    pivot_wider(names_from = Status, values_from = mean_value)

  # Model
  lm_model <- lm(Test ~ Train, data = lm_data)
  # Create threshold line
  thr = 1
  lm_data <- lm_data |>
    mutate(
      predicted = fitted(lm_model),  
      upper_threshold = predicted + thr,  
      lower_threshold = predicted - thr
    )

  # Correlation
  cor_pearson <- cor(lm_data$Train, lm_data$Test, method = "pearson")
  cor_spearman <- cor(lm_data$Train, lm_data$Test, method = "spearman")
  # P-value
  cor_p_pearson <- pvalue(independence_test(Test ~ Train, data = lm_data))
  cor_p_spearman <- pvalue(independence_test(rank(Test) ~ rank(Train), data = lm_data))
  # Dataframe
  test <- data.frame(
    Pearson = cor_pearson, 
    Spearman = cor_spearman,
    Prediction = "Subset",
    Status = "Test"
    )
  # Label
  pearson_label <- bquote(R[Pearson] == .(round(cor_pearson, 3)))
  spearman_label <- bquote(R[Spearman] == .(round(cor_spearman, 3)))

  # Train and test ancestry
  x_ancestry <- unique(filtered_data$Ancestry[filtered_data$Status == "Train"])
  y_ancestry <- unique(filtered_data$Ancestry[filtered_data$Status == "Test"])
  # Compared conditions
  condition_1 <- unique(raw_metric_dge$coef)[1]
  condition_2 <- unique(raw_metric_dge$coef)[2]

  # Plot
  prop_train_test <- lm_data |>
    ggplot(
      aes(
        x = Train,
        y = Test
      )
    ) +
    geom_point(
      size = point_size,
      color = "lightgrey"
    ) +
    geom_point(
      size = point_size,
      color = "lightgrey"
    ) +
    geom_point(
      aes(color = "dummy1"),
      size = NA
    ) +
    geom_point(
      aes(color = "dummy2"),
      size = NA
    ) +
    geom_smooth(
      method = "lm",
      color = "blue",
      linewidth = point_size,
      alpha = 0.5
    ) +
    geom_line(
      aes(x = Train, y = upper_threshold), 
      color = "blue", 
      linetype = "dashed",
      linewidth = point_size,
      alpha = 0.5
    ) + 
    geom_line(
      aes(x = Train, y = lower_threshold), 
      color = "blue", 
      linetype = "dashed",
      linewidth = point_size,
      alpha = 0.5
    ) +
    scale_color_manual(
      values = c(
        "dummy1" = "black",
        "dummy2" = "black"
        ),
      breaks = c(
        "dummy1",
        "dummy2"
        ),
      labels = c(
        pearson_label,
        spearman_label
      )
    ) +
    guides(
      color = guide_legend(
        title = "Subset",
        override.aes = list(size = 0, alpha = 0)
      )
    ) +
    labs(
      x = paste(x_ancestry, condition_1, "vs", condition_2, "(avg. logFC)"),
      y = paste(y_ancestry, condition_1, "vs", condition_2, "(avg. logFC)")
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend() +
    theme(
      legend.title = element_text(face = "bold"),
      legend.justification = c(0, 1),
      legend.position = c(0.05, 1.0), 
      legend.spacing.x = unit(0, "cm"), 
      legend.spacing.y = unit(0, "cm"),
      legend.background = element_rect(fill = "transparent"),  
      legend.key.spacing = unit(0, "cm")
    )

 
  # Bar plot for proportion
  small_bar_plot <- bind_rows(inference, test) |>
    pivot_longer(
      cols = c(
        Pearson,
        Spearman
      ),
      names_to = "Metric",
      values_to = "Value"
    ) |>
    ggplot(
      aes(
        x = fct_rev(Prediction),
        y = Value
      )
    ) +
    geom_bar(
      stat = "identity", 
      width = 0.7
    ) +
    facet_grid(
      rows = vars(Metric)
    ) +
    scale_x_discrete(
      labels = c(
        "Subset" = "S", 
        "Ancestry" = "A"
        )
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    ) +
    labs(
      x = "Prediction",
      y = "Y"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_white_strip() +
    theme_zero_margin() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),  
      panel.background = element_rect(fill = "transparent", color = NA),  
      axis.title = element_blank()
    )
  
  # Pacthwork for a tile
  insert_patchwork <- prop_train_test +
    inset_element(small_bar_plot, left = 0.7, bottom = 0.01, right = 1.1, top = 0.6) 

  tile_patchwork <- wrap_elements(
    (insert_patchwork + 
    prop_train_inference) +
    plot_annotation(title = inf_n_label) &
    theme(
        plot.title = element_text(hjust = 0.5, size = 5)
        )
    )

  # Add to list
  gg_proportion_list[[as.character(proportion)]] <- tile_patchwork
}
# Final layout (patchwork)
contrast_distribution_avg_logFC <- wrap_plots(gg_proportion_list, ncol = 2)
  
# Save
ggsave(filename = "Contrast_distribution_avg_logFC.pdf", 
      plot = contrast_distribution_avg_logFC, 
      path = path_to_save_location, 
      width = 8, height = 5
      )

# ---- Contrast density R-values (DGE) -----
