#' Filter Features by Variance and Sample Presence
#'
#' Filters features (e.g., genes or probes) based on variance across samples
#' and the number of samples with non-zero expression or signal.
#'
#' @param data A numeric matrix (samples x features), such as methylation beta values.
#' @param var_threshold Numeric. Minimum variance required to retain a feature. Default is 0.01.
#' @param min_samples_ratio Numeric. Minimum fraction (0–1) of samples where the feature must be non-zero. Default is 0.5.
#'
#' @return A logical vector indicating which features to retain.
#' @export
filter_by_variance <- function(
  data,
  var_threshold      = 0.01,
  min_samples_ratio  = 0.5
) {
  if (!is.matrix(data)) {
    stop("The 'data' argument must be a numeric matrix.")
  }

  if (!is.numeric(var_threshold) || length(var_threshold) != 1) {
    stop("The 'var_threshold' argument must be a single numeric value.")
  }

  if (!is.numeric(min_samples_ratio) || min_samples_ratio <= 0 || min_samples_ratio > 1) {
    stop("The 'min_samples_ratio' argument must be a numeric value between 0 and 1.")
  }

  num_samples   <- nrow(data)
  min_samples   <- floor(min_samples_ratio * num_samples)
  gene_variance <- apply(data, 2, stats::var, na.rm = TRUE)
  detected_samples <- colSums(data > 0, na.rm = TRUE)

  keep <- (gene_variance > var_threshold) &
          (detected_samples >= min_samples)

  return(keep)
}
