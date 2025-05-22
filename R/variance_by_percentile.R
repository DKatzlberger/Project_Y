#' Calculate Variance Threshold by Percentile
#'
#' Computes a threshold for feature-wise variance based on the given percentile.
#' Commonly used to filter features (e.g., genes, probes) with low variance.
#'
#' @param data A numeric matrix of expression or beta values (samples x features).
#' @param percentile Numeric. Percentile (0–100) to use for calculating the variance threshold.
#'        Default is 25.
#'
#' @return A single numeric value representing the variance threshold.
#' @export
variance_by_percentile <- function(
  data,
  percentile = 25
) {
  if (!is.matrix(data)) {
    stop("⚠ The 'data' argument must be a numeric matrix.")
  }

  if (!is.numeric(percentile) || length(percentile) != 1 || percentile < 0 || percentile > 100) {
    stop("⚠ The 'percentile' argument must be a number between 0 and 100.")
  }

  gene_variances <- apply(data, 2, stats::var, na.rm = TRUE)
  threshold <- stats::quantile(gene_variances, probs = percentile / 100)

  return(threshold)
}
