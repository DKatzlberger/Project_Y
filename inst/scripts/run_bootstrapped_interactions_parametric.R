#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────
# Cross-Ancestry Parametric Main Script (Development Mode)
# ───────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(anndata)
  library(yaml)
  library(glue)
  library(uuid)
  library(reticulate)
  library(parallel)
  library(purrr)
  library(furrr)
  library(limma)
  library(patchwork)
  library(ggplot2)

  use_condaenv(
    condaenv = "/opt/conda/envs/ancestry/bin/python",
    required = TRUE
  )
})

# Load package function files (Development Mode)
function_dir <- "R"
function_files <- list.files(
  path       = function_dir,
  full.names = TRUE,
  pattern    = "\\.R$"
)

suppressMessages(
  invisible(lapply(function_files, source))
)

if (length(function_files) == 0) {
  cat("⚠ No R function files found in 'R/' directory.\n")
} else {
  cat("✔ Loaded function files from R/: \n")
  cat(paste0("• ", basename(function_files), collapse = "\n"), "\n")
}


# YAML settings file (from CLI or fallback)
YAML_FILE <- yaml_argument_loader(
  default_yaml = "inst/extdata/example_settings_bootstrapped_interactions_parametric.yaml"
)

# Default settings file
DEFAULT_FILE <- get_default_settings_path(
  file_name = "default_settings_bootstrapped_interactions_parametric.yaml",
  base_dir  = file.path("inst", "config")
)

# Merge and load
setup <- load_and_merge_settings(
  default_path = DEFAULT_FILE,
  user_path    = YAML_FILE
)

# Check required settings
required_settings <- c(
  "output_column", 
  "class_0", 
  "class_1", 
  "ancestry_column", 
  "train_ancestry", 
  "infer_ancestry", 
  "data_path", 
  "tech", 
  "normalization_method",
  "output_directory"
)

check_required_settings(
  settings      = setup,
  required_keys = required_settings
)


# Create output directory and save finalized settings
path_to_save_location <- setup$output_directory
if (!dir.exists(path_to_save_location)) {
  dir.create(
    path      = path_to_save_location,
    recursive = TRUE
  )
}

save_settings(
  settings   = setup,
  file_path  = file.path(setup$output_directory, "Settings.yaml")
)

# Load and validate input data
data_path       <- setup$data_path
output_column   <- setup$output_column
ancestry_column <- setup$ancestry_column

is_h5ad_file(
  file_path = data_path
)

adata <- read_h5ad(
  filename = data_path
)

required_columns <- c(
  output_column,
  ancestry_column
)

check_columns_exist(
  df            = adata$obs,
  required_cols = required_columns
)


# Define and validate classification labels
class_0 <- setup$class_0
class_1 <- setup$class_1

comparison <- c(
  class_0,
  class_1
)

check_values_exist(
  df       = adata$obs,
  column   = output_column,
  expected = comparison
)

adata <- adata[
  adata$obs[[output_column]] %in% comparison,
]

adata$obs[[output_column]] <- factor(
  x      = adata$obs[[output_column]],
  levels = comparison
)

# Define and validate ancestry groups
train_ancestry  <- setup$train_ancestry
infer_ancestry  <- setup$infer_ancestry

ancestries <- c(
  train_ancestry,
  infer_ancestry
)

check_values_exist(
  df       = adata$obs,
  column   = ancestry_column,
  expected = ancestries
)

adata <- adata[
  adata$obs[[ancestry_column]] %in% ancestries,
]

adata$obs[[ancestry_column]] <- factor(
  x      = adata$obs[[ancestry_column]],
  levels = ancestries
)


# Sample size
p_count <- plot_output_column_count(
  data        = adata$obs,
  x           = ancestry_column,
  fill        = output_column,
  point_size  = 0.5
)

p_proportions <- plot_output_column_proportion(
  data        = adata$obs,
  x           = ancestry_column,
  fill        = output_column,
  point_size  = 0.5
)

p <- p_count + p_proportions + patchwork::plot_layout(guides = "collect") & 
  ggplot2::theme(legend.position = "bottom")

save_name <- file.path(path_to_save_location, "QC_sample_sizes.pdf")

size <- estimate_plot_size(p)

save_ggplot(
  plot      = p,
  save_path = save_name,
  width     = size$width,
  height    = size$height
)

# Correlation heatmap
save_name <- file.path(path_to_save_location, "QC_correaltion_heatmap.pdf")
plot_correlation_heatmap(
  expr_data       = t(adata$X),
  metadata        = adata$obs,
  ancestry_column = ancestry_column,
  output_column   = output_column,
  file_path       = save_name,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_names      = FALSE,
  width           = 5,
  height          = 5
)


# Feature filtering
filter_features  <- setup$filter_features
tech             <- setup$tech
percentile       <- setup$percentile
data_type        <- setup$data_type

save_name <- file.path(path_to_save_location, "QC_mean_variance_trend.pdf")
if (filter_features && tech == "transcriptomics") {
  filtered_data <- filter_features_transcriptomics(
    adata      = adata,
    percentile = percentile,
    data_type  = data_type,
    plot_path  = save_name
  )
} else{
  stop("Filtering did not work!")
}

# Normalization
filtered_data <- voom_normalization(
  adata = filtered_data
)

p_before <- plot_density_of_counts(
  filtered_data$raw$X, 
  x_axis_label = setup$data_type
  )
p_after <- plot_density_of_counts(
  filtered_data$voom$X, 
  x_axis_label = "logCPM"
  )

p <- p_before + p_after + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save
save_name <- file.path(path_to_save_location, "QC_density_counts.pdf")
size <- estimate_plot_size(p)
save_ggplot(
  plot = p,
  save_path = save_name,
  width = size$width,
  height = size$height
)

# Collect results from many subsets
all_interactions <- list()

# Seed loop over randomly chosen seeds
seed <- setup$reproducibility_key
n_subsets <- setup$n_subsets
generate_random_integers <- function(seed, n = 10, min = 1, max = 1e6) {
  # Save current RNG state
  old_rng <- .Random.seed
  rng_kind <- RNGkind()
  on.exit({ RNGkind(rng_kind[1], rng_kind[2], rng_kind[3]); .Random.seed <<- old_rng })
  # Use local seed
  set.seed(seed)
  sample(min:max, n, replace = FALSE)
}

random_seed <- generate_random_integers(seed = seed, n = n_subsets)

for (seed in random_seed){

  # Create a subfolder seed
  seed_dir <- file.path(path_to_save_location, paste0("subset_", seed))
  if (!dir.exists(seed_dir)) dir.create(seed_dir, recursive = TRUE)

  # Ancestry stratification
  save_name <- file.path(seed_dir, "QC_ancestry_stratification.pdf")
  strat_result <- stratify_ancestry_subsets(
    adata = filtered_data$voom,
    ancestry_column = ancestry_column,
    output_column = output_column,
    train_ancestry = train_ancestry,
    infer_ancestry = infer_ancestry,
    seed = seed,
    plot_path = save_name
    )

  test_adata  <- strat_result$test_adata
  infer_adata <- strat_result$infer_adata
  train_adata <- strat_result$train_adata

  # Bootstrap interaction
  interaction_res <- bootstrap_interaction(
    X = test_adata$X,            
    Y = infer_adata$X,            
    MX = test_adata$obs,          
    MY = infer_adata$obs,         
    g_col = output_column,        
    a_col = ancestry_column,      
    B = setup$n_bootstraps,                    
    seed = seed                   
    )

  # Save for across seeds
  bootstrapped_interaction <- interaction_res$observed
  all_interactions[[as.character(seed)]] <- bootstrapped_interaction
  # Save interaction statistic
  write.csv(
    bootstrapped_interaction, 
    file = file.path(seed_dir, "Bootstrapped_interaction.csv")
    )

  # Volcano plot of bootstrapped interaction
  p <- plot_volcano(
    data = bootstrapped_interaction,
    effect_col = "T_obs",
    p_value_col = "adj_p_value",
    effect_label = "logFC",
    title = paste("Bootstrapped:", unique(bootstrapped_interaction$coef))
    )

  # Save volcano plot of bootstrapped interaction
  save_name <- file.path(seed_dir, "Bootstrapped_volcano_plot.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot = p,
    save_path = save_name,
    width = size$width,
    height = size$height
    )

  # Get significant interaction
  sig_bootstrapped <- bootstrapped_interaction[
    abs(T_obs) > 1 & adj_p_value < 0.05, feature]

  if (length(sig_bootstrapped) != 0){
    # Prepare z-score values
    train_zscore <- prepare_normalized_data(
      matrix = t(train_adata$X),
      meta = strat_result$train_adata$obs,
      features = sig_bootstrapped
      )[, Set := "Train"]

    test_zscore <- prepare_normalized_data(
      matrix = t(test_adata$X),
      meta = strat_result$test_adata$obs,
      features = sig_bootstrapped
      )[, Set := "Test"]

    infer_zscore <- prepare_normalized_data(
      matrix = t(infer_adata$X),
      meta = strat_result$infer_adata$obs,
      features = sig_bootstrapped
      )[, Set := "Infer"]

    # Prepare factor for different sets
    zscore <- rbind(train_zscore, test_zscore, infer_zscore)
    zscore[, Set := factor(Set, levels = c("Train", "Test", "Infer"))]

    # Boxplot of expression zscores
    p <- plot_expression_boxplots(
      data = zscore,
      x = "Set",
      y = "zscore",
      fill = output_column
      )

    # Save the boxplots
    save_name <- file.path(seed_dir, "Bootstrapped_sig_features.pdf")
    size <- estimate_plot_size(p)
    save_ggplot(
      plot = p,
      save_path = save_name,
      width = size$width + 2,
      height = size$height + 2
      )
  }

  # Limma quality control
  limma_res <- do_limma_interactions(
    Y = rbind(test_adata$X, infer_adata$X),           
    M = rbind(test_adata$obs, infer_adata$obs),   
    g_col = output_column,          
    a_col = ancestry_column,            
    a_train = train_ancestry,     
    a_infer = infer_ancestry,     
    g_levels = comparison,   
    normalize = FALSE
    )

  # Get top table of interactions
  limma_interaction <- extract_top_table(
    limma_fit = limma_res$fit_contrasts,
    feature_column = "feature",
    exclude_intercept = TRUE
  )[coef == "deltaB_eur - deltaB_admix"]

  # Volcano plot of interactions
  p <- plot_volcano(
    data = limma_interaction,
    effect_col = "logFC",
    p_value_col = "adj.P.Val",
    effect_label = "logFC",
    title = paste("Limma:", unique(limma_interaction$coef))
  )

  # Save volcano plot of limma interaction
  save_name <- file.path(seed_dir, "Limma_volcano_plot.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot = p,
    save_path = save_name,
    width = size$width,
    height = size$height
    )

  # Validate pvalues with limma
  validate_pvals <- merge(
    x = limma_interaction[, .(feature, p_limma = P.Value, adj_p_limma = adj.P.Val)],
    y = bootstrapped_interaction[, .(feature, p_boot = p_value, adj_p_boot = adj_p_value)],
    by = "feature"
    )

  # Transform to be more linear
  validate_pvals[, `:=`(
    log_p_limma = -log10(p_limma),
    log_p_boot = -log10(p_boot))
    ]

  # Sigificance source
  validate_pvals[, signif_source := fifelse(
    adj_p_limma < 0.05 & adj_p_boot < 0.05, "Both",
    fifelse(adj_p_limma < 0.05, "Limma only",
    fifelse(adj_p_boot  < 0.05, "Bootstrap only", "None")))
    ]

  # Correlation of pvalues (validation)
  p <- plot_pvalue_concordance(
    data   = validate_pvals,
    p_col_x = "p_limma",
    p_col_y = "p_boot",
    signif_col = "signif_source",
    method_x_label = "limma",
    method_y_label = "bootstrap",
    log_cap = 5,
    title = "Limma vs Bootstrap: Pvalue concordance"
    )

  # Save validation of pvalues 
  save_name <- file.path(seed_dir, "Pvalue_concordance.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot = p,
    save_path = save_name,
    width = size$width,
    height = size$height
    )

  # Bootstrap pearson correlation
  pearson_res <- bootstrap_correlation(
    X = test_adata$X,             
    Y = infer_adata$X,           
    Z = train_adata$X,            
    MX = test_adata$obs,          
    MY = infer_adata$obs,         
    MZ = train_adata$obs,         
    g_col = output_column,        
    a_col = ancestry_column,      
    method = "pearson",           
    B = setup$n_bootstraps,                     
    seed = seed                   
    )
  
  # Save pearson statistic
  bootstrapped_pearson <- pearson_res$observed
  write.csv(
    bootstrapped_pearson, 
    file = file.path(seed_dir, "Bootstrapped_pearson.csv")
    )

  # Bootstrap spearman correlation
  spearman_res <- bootstrap_correlation(
    X = test_adata$X,             
    Y = infer_adata$X,           
    Z = train_adata$X,            
    MX = test_adata$obs,          
    MY = infer_adata$obs,         
    MZ = train_adata$obs,         
    g_col = output_column,        
    a_col = ancestry_column,      
    method = "spearman",           
    B = setup$n_bootstraps,                     
    seed = seed                   
    )
  
  # Save spearman statistics
  bootstrapped_spearman <- spearman_res$observed
  write.csv(
    bootstrapped_spearman, 
    file = file.path(seed_dir, "Bootstrapped_spearman.csv")
    )

  # Prepare raw correlation 
  pearson_res$correlation$p_value <- bootstrapped_pearson$p_value
  spearman_res$correlation$p_value <- bootstrapped_spearman$p_value
  
  p <- plot_correlation(
    data = rbind(pearson_res$correlation, spearman_res$correlation),
    x_col = "target",
    y_col = "correlation",
    facet_col = "method",
    pval_col = "p_value",
    ci_lower_col = "ci_lower",
    ci_upper_col = "ci_upper",
    point_size = 0.5,
    title = NULL
    )

  # Save correlation bar plot
  save_name <- file.path(seed_dir, "Boostrapped_correlation_plot.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot = p,
    save_path = save_name,
    width = size$width,
    height = size$height
    )
}

# Summaries across seeds
all_interactions_DT <- rbindlist(all_interactions, idcol = "seed")
summary_interactions <- all_interactions_DT[, .(
  mean_T_obs = mean(T_obs, na.rm = TRUE),
  sd_T_obs = sd(T_obs, na.rm = TRUE),
  median_p = median(p_value, na.rm = TRUE),
  prop_signif = mean(p_value < 0.05, na.rm = TRUE),
  n_seeds = .N
), by = feature]

# Adjustment of pvalues
summary_interactions[, adj_median_p := p.adjust(median_p, method = "BH")]

# Save the summarized statistics
fwrite(
  summary_interactions, 
  file = file.path(path_to_save_location, "Bootstrapped_interaction.csv")
  )

# Visualize median pvalue with reproducibility 
p <- plot_volcano(
  data = summary_interactions,
  effect_col = "mean_T_obs",
  p_value_col = "median_p",
  effect_label = "Mean logFC",
  color_by = "prop_signif",
  title = "Volcano Plot: Mean effect size vs. reproducibility"
  )

# Save correlation bar plot
save_name <- file.path(path_to_save_location, "Weighted_volcano_plot.pdf")
size <- estimate_plot_size(p)
save_ggplot(
  plot = p,
  save_path = save_name,
  width = size$width,
  height = size$height
  )

# Top hits with this method
top_hits <- summary_interactions[
  median_p < 0.05 &
  prop_signif >= 0.8 &
  abs(mean_T_obs) > 1,
  head(feature, 10)]

# Prepare expression zscore
t_a <- filtered_data$voom[filtered_data$voom$obs[[ancestry_column]] %in% train_ancestry]
train_zscore <- prepare_normalized_data(
  matrix = t(t_a$X),
  meta = t_a$obs,
  features = top_hits
  )[, Set := train_ancestry]

i_a <- filtered_data$voom[filtered_data$voom$obs[[ancestry_column]] %in% infer_ancestry]
infer_zscore <- prepare_normalized_data(
  matrix = t(i_a$X),
  meta = i_a$obs,
  features = top_hits
  )[, Set := infer_ancestry]

# Prepare factor for different sets
zscore <- rbind(train_zscore, infer_zscore)
zscore[, Set := factor(Set, levels = c(train_ancestry, infer_ancestry))]

# Boxplot of expression zscores
p <- plot_expression_boxplots(
  data = zscore,
  x = "Set",
  y = "zscore",
  fill = output_column
  )

# Save the boxplots
save_name <- file.path(path_to_save_location, "Bootstrapped_top_hits.pdf")
size <- estimate_plot_size(p)
save_ggplot(
  plot = p,
  save_path = save_name,
  width = size$width + 2,
  height = size$height + 2
  )


# Limma on all data
limma_res <- do_limma_interactions(
    Y = filtered_data$voom$X,           
    M = filtered_data$voom$obs,   
    g_col = output_column,          
    a_col = ancestry_column,            
    a_train = train_ancestry,     
    a_infer = infer_ancestry,     
    g_levels = comparison,   
    normalize = FALSE
    )

# Get top table of interactions
limma_interaction <- extract_top_table(
  limma_fit = limma_res$fit_contrasts,
  feature_column = "feature",
  exclude_intercept = TRUE
  )[coef == "deltaB_eur - deltaB_admix"]

# Volcano plot of interactions
p <- plot_volcano(
  data = limma_interaction,
  effect_col = "logFC",
  p_value_col = "adj.P.Val",
  effect_label = "logFC",
  title = paste("Limma:", unique(limma_interaction$coef))
  )

# Save volcano plot of limma interaction
save_name <- file.path(path_to_save_location, "Limma_volcano_plot.pdf")
size <- estimate_plot_size(p)
save_ggplot(
  plot = p,
  save_path = save_name,
  width = size$width,
  height = size$height
  )

# Top hits limma
top_hits <- limma_interaction[
  adj.P.Val < 0.05 & 
  abs(logFC) > 1,
  head(feature, 10)]

# Prepare expression zscore
t_a <- filtered_data$voom[filtered_data$voom$obs[[ancestry_column]] %in% train_ancestry]
train_zscore <- prepare_normalized_data(
    matrix = t(t_a$X),
    meta = t_a$obs,
    features = top_hits
    )[, Set := train_ancestry]

i_a <- filtered_data$voom[filtered_data$voom$obs[[ancestry_column]] %in% infer_ancestry]
infer_zscore <- prepare_normalized_data(
  matrix = t(i_a$X),
  meta = i_a$obs,
  features = top_hits
  )[, Set := infer_ancestry]

# Prepare factor for different sets
zscore <- rbind(train_zscore, infer_zscore)
zscore[, Set := factor(Set, levels = c(train_ancestry, infer_ancestry))]

# Boxplot of expression zscores
p <- plot_expression_boxplots(
  data = zscore,
  x = "Set",
  y = "zscore",
  fill = output_column
  )

# Save the boxplots
save_name <- file.path(path_to_save_location, "Limma_top_hits.pdf")
size <- estimate_plot_size(p)
save_ggplot(
  plot = p,
  save_path = save_name,
  width = size$width + 2,
  height = size$height + 2
  )
