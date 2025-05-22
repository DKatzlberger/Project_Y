#' Filter Transcriptomic Features by Signal and Variance
#'
#' Filters transcriptomic features (e.g., RNA-seq counts) based on signal
#' (e.g., logCPM) and variance thresholds defined by percentiles. If a
#' `plot_path` is provided, a before/after mean-variance plot is saved.
#'
#' @param adata An AnnData-like object with a matrix `adata$X` (samples x features).
#' @param percentile Numeric. Percentile threshold (0–100) for filtering signal and variance.
#' @param plot_path Optional. File path to save the before/after QC plot.
#' @param data_type Optional. Label for the x-axis in the plot. Required if `plot_path` is given.
#'
#' @return A filtered AnnData-like object with low-signal and low-variance features removed.
#' @export
filter_features_transcriptomics <- function(
  adata,
  percentile,
  plot_path = NULL,
  data_type = NULL
) {
  # Step 1: Filter by signal
  norm_factors <- calculate_tmm_norm_factors(adata$X)
  log_cpm      <- cpm(adata$X, norm_factors = norm_factors, log = TRUE)

  min_signal   <- signal_by_percentile(log_cpm, percentile)
  keep_signal  <- filter_by_signal(log_cpm, min_signal)
  filtered     <- adata[, keep_signal]

  # Step 2: Filter by variance
  norm_factors <- calculate_tmm_norm_factors(filtered$X)
  log_cpm      <- cpm(filtered$X, norm_factors = norm_factors, log = TRUE)

  min_variance <- variance_by_percentile(log_cpm, percentile)
  keep_var     <- filter_by_variance(log_cpm, var_threshold = min_variance)
  filtered     <- filtered[, keep_var]

  # Step 3: Optional plot
  if (!is.null(plot_path)) {
    if (is.null(data_type)) {
      stop("If 'plot_path' is provided, 'data_type' must also be specified.")
    }

    data_before <- log2(adata$X + 0.5)
    data_after  <- log2(filtered$X + 0.5)
    x_axis      <- paste0("log2(", data_type, " + 0.5)")

    p_before <- plot_mean_variance_trend(data_before, x_axis_label = x_axis) +
      ggplot2::ggtitle("Before filtering")
    p_after <- plot_mean_variance_trend(data_after, x_axis_label = x_axis) +
      ggplot2::ggtitle("After filtering")

    plot <- p_before / p_after
    size <- estimate_plot_size(plot)

    save_ggplot(
      plot      = plot,
      save_path = plot_path,
      width     = size$width,
      height    = size$height
    )
  }

  return(filtered)
}
