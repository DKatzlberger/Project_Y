#' Run the cross-ancestry parametric pipeline
#'
#' @return No return value. Side effects include saving plots and result files.
#' @export
bootstrapped_interactions_parametric_pipeline <- function(
  seed,
  output_column,
  class_0,
  class_1,
  ancestry_column,
  train_ancestry,
  infer_ancestry,
  data_path,
  tech,
  normalization_method,
  output_directory,
  filter_features,
  percentile,
  data_type,
  n_bootstraps,
  vis_features,
  sanity_check,
  ...
) {
  # Load necessary libraries and configure conda environment
  suppressPackageStartupMessages({
    library(data.table)
    library(anndata)
    library(glue)
    library(uuid)
    library(reticulate)
    library(parallel)
    library(purrr)
    library(furrr)
    library(limma)
    library(patchwork)

    use_condaenv(
      condaenv = "/opt/conda/envs/ancestry/bin/python",
      required = TRUE
    )
  })

  source("R/voom_normalization.R")

  # Create output directory and save settings if required
  path_to_save_location <- output_directory
    if (!dir.exists(path_to_save_location)) {
      dir.create(path_to_save_location, recursive = TRUE)
    }

  setup <- list(
      seed                 = seed,
      output_column        = output_column,
      class_0              = class_0,
      class_1              = class_1,
      ancestry_column      = ancestry_column,
      train_ancestry       = train_ancestry,
      infer_ancestry       = infer_ancestry,
      data_path            = data_path,
      tech                 = tech,
      normalization_method = normalization_method,
      output_directory     = output_directory,
      filter_features      = filter_features,
      percentile           = percentile,
      data_type            = data_type,
      n_bootstraps         = n_bootstraps,
      vis_features         = vis_features,
      sanity_check         = sanity_check
    )

    save_settings(
      setup,
      file.path(path_to_save_location, "Settings.yaml")
    )


  # Load input data and filter to relevant samples
  is_h5ad_file(data_path)
  adata <- read_h5ad(data_path)

  check_columns_exist(
    df            = adata$obs, 
    required_cols = c(output_column, ancestry_column)
    )

  comparison <- c(class_0, class_1)
  check_values_exist(
    df       = adata$obs, 
    column   = output_column, 
    expected = comparison
    )

  adata <- adata[adata$obs[[output_column]] %in% comparison, ]
  adata$obs[[output_column]] <- factor(
    x      = adata$obs[[output_column]], 
    levels = comparison
    )

  ancestries <- c(train_ancestry, infer_ancestry)
  check_values_exist(
    df       = adata$obs, 
    column   = ancestry_column, 
    expected = ancestries
    )

  adata <- adata[adata$obs[[ancestry_column]] %in% ancestries, ]
  adata$obs[[ancestry_column]] <- factor(
    x      = adata$obs[[ancestry_column]], 
    levels = ancestries
    )

  # Set seed for reproducibility
  set.seed(seed)

  # Plot sample sizes and proportions
  p_count <- plot_output_column_count(
    data       = adata$obs, 
    x          = ancestry_column, 
    fill       = output_column, 
    point_size = 0.5
    )

  p_proportions <- plot_output_column_proportion(
    data       = adata$obs, 
    x          = ancestry_column, 
    fill       = output_column, 
    point_size = 0.5
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

  # Plot correlation heatpmap
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


  # Filter features if specified
  if (filter_features && tech == "transcriptomics") {
    filtered_data <- filter_features_transcriptomics(
      adata      = adata,
      percentile = percentile,
      data_type  = data_type,
      plot_path  = file.path(path_to_save_location, "QC_mean_variance_trend.pdf") 
    )

  } else {
    stop("Filtering did not work!")
  }
  # Save features 
  save_name <- file.path(path_to_save_location, "Features.yaml")
  write_yaml(filtered_data$var_names, save_name)

  setup$n_features <- length(filtered_data$var_names)
  save_settings(
    settings   = setup,
    file_path  = file.path(output_directory, "Settings.yaml")
  )

  # Stratify data by ancestry for training and inference
  strat_result <- stratify_ancestry_subsets(
    adata           = filtered_data,
    ancestry_column = ancestry_column,
    output_column   = output_column,
    train_ancestry  = train_ancestry,
    infer_ancestry  = infer_ancestry,
    seed            = seed,
    plot_path       = file.path(path_to_save_location, "QC_ancestry_stratification.pdf")
  )

  # Calculate observed logFC values
  normalization_fn <- match.fun(normalization_method)
  train_obs <- calculate_logfc(
    matrix        = strat_result$train_adata$X, 
    meta          = strat_result$train_adata$obs, 
    group_column  = output_column, 
    normalization = normalization_fn
    )

  test_obs  <- calculate_logfc(
    matrix        = strat_result$test_adata$X, 
    meta          = strat_result$test_adata$obs, 
    group_column  = output_column, 
    normalization = normalization_fn
    )

  infer_obs <- calculate_logfc(
    matrix        = strat_result$infer_adata$X, 
    meta          = strat_result$infer_adata$obs, 
    group_column  = output_column, 
    normalization = normalization_fn
    )

  # Compute observed statistics for differences and correlations
  interaction_obs <- compute_logfc_diff_obs(
    test    = test_obs$logfc, 
    infer   = infer_obs$logfc,
    join_by = "feature"
    )

  pearson_obs <- compute_logfc_correlation_difference(
    train  = train_obs$logfc, 
    test   = test_obs$logfc, 
    infer  = infer_obs$logfc, 
    method = "pearson"
    )

  spearman_obs <- compute_logfc_correlation_difference(
    train  = train_obs$logfc, 
    test   = test_obs$logfc, 
    infer  = infer_obs$logfc, 
    method = "spearman"
    )



  # Bootstrap samples from null hypothesis distribution
  H0_matrix  <- rbind(strat_result$test_adata$X, strat_result$infer_adata$X)
  H0_meta    <- rbind(strat_result$test_adata$obs, strat_result$infer_adata$obs)

  test_boot <- bootstrap_logfc(
    matrix        = H0_matrix, 
    meta          = H0_meta, 
    size          = nrow(strat_result$test_adata), 
    group_column  = output_column, 
    normalization = normalization_fn, 
    n_iterations  = n_bootstraps
    )

  infer_boot <- bootstrap_logfc(
    matrix        = H0_matrix, 
    meta          = H0_meta, 
    size          = nrow(strat_result$infer_adata), 
    group_column  = output_column, 
    normalization = normalization_fn, 
    n_iterations  = n_bootstraps
    )

  # Compute bootstrap statistics
  interaction_boot <- compute_logfc_diff_boot(
    test     = test_boot, 
    infer    = infer_boot,
    join_by  = "feature",
    boot_col = "bootstrap"
    )

  train_boot <- train_obs$logfc[, .(bootstrap = unique(test_boot$bootstrap)), by = .(feature, logFC)]

  pearson_boot <- compute_logfc_correlation_difference(
    train  = train_boot, 
    test   = test_boot, 
    infer  = infer_boot, 
    method = "pearson"
    )

  spearman_boot <- compute_logfc_correlation_difference(
    train  = train_boot, 
    test   = test_boot, 
    infer  = infer_boot, 
    method = "spearman"
    )



  # Define the features to visualize
  if (!is.null(vis_features)){
    vis_features <- vis_features
  } else{
    vis_features <- interaction_obs[1:9, feature]
  }


  # Calculate p-values and adjust
  interaction_summary <- calculate_pvalues(
    observed_dt  = interaction_obs,
    bootstrap_dt = interaction_boot,
    value_col    = "logFC",
    by_cols      = "feature",
    is_global    = FALSE,
    plot_path    = file.path(path_to_save_location, "QC_null_interaction.pdf"),
    features     = vis_features
  )
  interaction_summary[, p_adj_emp := p.adjust(p_emp, method = "BH")]
  interaction_summary[, p_adj_param := p.adjust(p_param, method = "BH")]


  # Pearson summary
  pearson_summary <- calculate_pvalues(
    observed_dt  = pearson_obs,
    bootstrap_dt = pearson_boot,
    value_col    = "difference",
    by_cols      = "feature",
    is_global    = TRUE,
    plot_path    = file.path(path_to_save_location, "QC_null_pearson.pdf")
  )


  # Spearman summary
  spearman_summary <- calculate_pvalues(
    observed_dt  = spearman_obs,
    bootstrap_dt = spearman_boot,
    value_col    = "difference",
    by_cols      = "feature",
    is_global    = TRUE,
    plot_path    = file.path(path_to_save_location, "QC_null_spearman.pdf") 
  )



  # Save summaries
  fwrite(interaction_summary, file.path(path_to_save_location, "Interaction_summary.csv"))

  correlation_summary <- rbind(pearson_summary, spearman_summary)
  fwrite(correlation_summary, file.path(path_to_save_location, "Correlation_summary.csv"))

  # Save significant features
  sig_features <- interaction_summary[abs(logFC) > 1 & p_adj_param < 0.05]$feature
  save_name    <- file.path(path_to_save_location, "Significant_features.yaml")
  write_yaml(sig_features, save_name)




  # Visualize the results
  path_to_save_location <- file.path(output_directory, "Visual_val")
  if (!dir.exists(path_to_save_location)) {
    dir.create(
      path      = path_to_save_location,
      recursive = TRUE
    )
  }

  # Interactions
  p <- plot_volcano(
    data         = interaction_summary,
    effect_col   = "logFC",
    p_value_col  = "p_adj_param",
    effect_label = "logFC",
    title        = NULL
  )

  # Save volcano plot
  save_name <- file.path(path_to_save_location, "Volcano_plot.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )



  # Correlations 
  p <- plot_correlation(
    data        = correlation_summary,
    test_label  = "Test",
    infer_label = "Inference",
    title       = NULL
    )
  
  # Save correlation bar plot
  save_name <- file.path(path_to_save_location, "Correlation_plot.pdf")
  size <- estimate_plot_size(p)
  save_ggplot(
    plot      = p,
    save_path = save_name,
    width     = size$width,
    height    = size$height
  )



  if ((!is.null(vis_features) && length(vis_features) > 0) ||
    (length(sig_features) > 0)) {

    # Determine which features to visualize
    if (!is.null(vis_features) && length(sig_features) > 0) {
      vis_features <- union(vis_features, sig_features)[1:9]
    } else if (length(sig_features) > 0) {
      vis_features <- sig_features[1:9]
    } else {
      vis_features <- vis_features[1:9]
    }

    # Prepare z-score data
    train_zscore <- prepare_normalized_data(
      matrix   = t(train_obs$norm_matrix),
      meta     = strat_result$train_adata$obs,
      features = vis_features
    )
    train_zscore[, Set := "Train"]

    test_zscore <- prepare_normalized_data(
      matrix   = t(test_obs$norm_matrix),
      meta     = strat_result$test_adata$obs,
      features = vis_features
    )
    test_zscore[, Set := "Test"]

    infer_zscore <- prepare_normalized_data(
      matrix   = t(infer_obs$norm_matrix),
      meta     = strat_result$infer_adata$obs,
      features = vis_features
    )
    infer_zscore[, Set := "Infer"]

    zscore <- rbind(train_zscore, test_zscore, infer_zscore)

    p <- plot_expression_boxplots(
      data = zscore,
      x    = "Set",
      y    = "zscore",
      fill = output_column
    )

    save_name <- file.path(path_to_save_location, "Interaction_boxplot.pdf")
    size      <- estimate_plot_size(p)
    save_ggplot(
      plot      = p,
      save_path = save_name,
      width     = size$width,
      height    = size$height
    )
  }


  # Optional: validate using limma interaction model
  if (sanity_check) {

    path_to_save_location <- file.path(output_directory, "Sanity_check_limma")
    if (!dir.exists(path_to_save_location)) {
      dir.create(
        path      = path_to_save_location,
        recursive = TRUE
      )
    }

    limma_matrix <- rbind(strat_result$test_adata$X, strat_result$infer_adata$X)
    limma_meta   <- rbind(strat_result$test_adata$obs, strat_result$infer_adata$obs)

    limma_res <- do_limma_interactions(
      expr_data       = t(limma_matrix),
      sample_meta     = limma_meta,
      ancestry_column = ancestry_column,
      output_column   = output_column,
      train_ancestry  = train_ancestry,
      infer_ancestry  = infer_ancestry,
      comparison      = comparison
    )

    fwrite(extract_top_table(limma_res$fit_means, "feature"), file.path(path_to_save_location, "Limma_means.csv"))
    fwrite(extract_top_table(limma_res$fit_contrasts, "feature"), file.path(path_to_save_location, "Limma_contrast.csv"))

    limma_contrasts <- extract_top_table(
      limma_fit      = limma_res$fit_contrasts,
      feature_column = "feature"
    )

    limma_interactions <- limma_contrasts[grepl("Interaction", coef)]

    p <- plot_volcano(
      data         = limma_interactions,
      effect_col   = "logFC",
      p_value_col  = "adj.P.Val",
      effect_label = "logFC",
      title        = NULL
    )

    save_name <- file.path(path_to_save_location, "Volcano_plot.pdf")

    size <- estimate_plot_size(p)

    save_ggplot(
      plot      = p,
      save_path = save_name,
      width     = size$width,
      height    = size$height
    )

  }
}
