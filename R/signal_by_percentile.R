#' Calculate Signal Threshold by Percentile
#'
#' Computes a signal threshold based on the sum of expression values per feature
#' (e.g., gene), using a specified percentile. Useful for filtering low-expression features.
#'
#' @param data A numeric matrix of expression values (samples x features).
#' @param percentile Numeric. Percentile (0–100) used to compute the signal threshold.
#'
#' @return A single numeric value representing the signal threshold.
#' @export
signal_by_percentile <- function(
  data,
  percentile
) {
  if (!is.numeric(percentile) || length(percentile) != 1 || percentile < 0 || percentile > 100) {
    stop("⚠ The `percentile` must be a number between 0 and 100.")
  }

  if (!is.matrix(data)) {
    stop("⚠ The `data` must be a numeric matrix.")
  }

  # Compute total signal for each feature (column sums)
  gene_signal <- colSums(data, na.rm = TRUE)

  # Return threshold at given percentile
  threshold <- stats::quantile(gene_signal, probs = percentile / 100)

  return(threshold)
}
