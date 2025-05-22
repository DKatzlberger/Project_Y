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
  filter_features   = TRUE,
  percentile        = 25,
  data_type         = NULL,
  n_bootstraps      = 1000,
  sanity_check      = TRUE,
  save_outputs      = TRUE
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
  if (isTRUE(save_outputs)) {
    if (!dir.exists(path_to_save_location)) {
      dir.create(path_to_save_location, recursive = TRUE)
    }
    save_settings(
      list(
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
        sanity_check         = sanity_check
      ),
      file.path(path_to_save_location, "Settings.yaml")
    )
  }

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

  # Plot sample sizes and proportions if required
  if (isTRUE(save_outputs)) {

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

    size <- estimate_plot_size(p)

    save_ggplot(p, file.path(path_to_save_location, "QC_sample_sizes.pdf"), size$width, size$height)

  }

  # Filter features if specified
  if (filter_features && tech == "transcriptomics") {
    filtered_data <- filter_features_transcriptomics(
      adata      = adata,
      percentile = percentile,
      data_type  = data_type,
      plot_path  = if (isTRUE(save_outputs)) file.path(path_to_save_location, "QC_mean_variance_trend.pdf") else NULL
    )

  } else {
    stop("Filtering did not work!")
  }

  # Stratify data by ancestry for training and inference
  strat_result <- stratify_ancestry_subsets(
    adata           = filtered_data,
    ancestry_column = ancestry_column,
    output_column   = output_column,
    train_ancestry  = train_ancestry,
    infer_ancestry  = infer_ancestry,
    seed            = seed,
    plot_path       = if (isTRUE(save_outputs)) file.path(path_to_save_location, "QC_ancestry_stratification.pdf") else NULL
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
  interaction_obs <- compute_logfc_difference(
    test  = test_obs, 
    infer = infer_obs
    )

  pearson_obs <- compute_logfc_correlation_difference(
    train  = train_obs, 
    test   = test_obs, 
    infer  = infer_obs, 
    method = "pearson"
    )

  spearman_obs <- compute_logfc_correlation_difference(
    train  = train_obs, 
    test   = test_obs, 
    infer  = infer_obs, 
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
  interaction_boot <- compute_logfc_difference(
    test  = test_boot, 
    infer = infer_boot
    )

  train_boot <- train_obs[, .(bootstrap = unique(test_boot$bootstrap)), by = .(feature, logFC)]

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

  # Calculate p-values and adjust
  interaction_summary <- calculate_pvalues(
    observed_dt  = interaction_obs,
    bootstrap_dt = interaction_boot,
    value_col    = "logFC",
    by_cols      = "feature",
    is_global    = FALSE,
    plot_path    = if (isTRUE(save_outputs)) file.path(path_to_save_location, "QC_null_interaction.pdf") else NULL,
    features     = interaction_obs[1:9, feature]
  )
  interaction_summary[, p_adj_emp := p.adjust(p_emp, method = "BH")]
  interaction_summary[, p_adj_param := p.adjust(p_param, method = "BH")]

  pearson_summary <- calculate_pvalues(
    observed_dt  = pearson_obs,
    bootstrap_dt = pearson_boot,
    value_col    = "difference",
    by_cols      = "feature",
    is_global    = TRUE,
    plot_path    = if (isTRUE(save_outputs)) file.path(path_to_save_location, "QC_null_pearson.pdf") else NULL
  )

  spearman_summary <- calculate_pvalues(
    observed_dt  = spearman_obs,
    bootstrap_dt = spearman_boot,
    value_col    = "difference",
    by_cols      = "feature",
    is_global    = TRUE,
    plot_path    = if (isTRUE(save_outputs)) file.path(path_to_save_location, "QC_null_spearman.pdf") else NULL
  )

  # Save result summaries if required
  if (isTRUE(save_outputs)) {
    fwrite(interaction_summary, file.path(path_to_save_location, "Interaction_summary.csv"))

    correlation_summary <- rbind(pearson_summary, spearman_summary)
    fwrite(correlation_summary, file.path(path_to_save_location, "Correlation_summary.csv"))
  }

  # Optional: validate using limma interaction model
  if (sanity_check && isTRUE(save_outputs)) {
    limma_matrix <- rbind(strat_result$test_adata$X, strat_result$infer_adata$X)
    limma_meta   <- rbind(strat_result$test_adata$obs, strat_result$infer_adata$obs)

    limma_res <- do_interactions(
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
  }
}
