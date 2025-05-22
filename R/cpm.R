#' Counts Per Million (CPM) Transformation
#'
#' Converts a matrix of raw counts to CPM (Counts Per Million) values.
#' Optionally applies normalization factors and a log2 transformation.
#'
#' @param data A numeric matrix of raw counts (samples x features).
#' @param norm_factors Optional numeric vector of normalization factors for each sample.
#'        Must match the number of rows in `data`.
#' @param log Logical. Whether to apply log2(CPM + 1) transformation. Default is TRUE.
#'
#' @return A numeric matrix of CPM or log-CPM values (same dimensions as `data`).
#' @export
cpm <- function(
  data,
  norm_factors = NULL,
  log          = TRUE
) {
  # Apply normalization factors if provided
  if (!is.null(norm_factors)) {
    if (length(norm_factors) != nrow(data)) {
      stop("⚠ Length of `norm_factors` must match number of samples (rows in data).")
    }

    if (any(norm_factors == 0)) {
      stop("⚠ Normalization factors contain zeros, which would cause division by zero.")
    }

    data <- sweep(data, 1, norm_factors, FUN = "/")
  }

  # Check for non-finite values
  if (any(!is.finite(data))) {
    stop("⚠ Input data contains NA, NaN, or Inf values.")
  }

  # Calculate total reads per sample
  total_counts <- rowSums(data, na.rm = TRUE)

  if (any(is.na(total_counts))) {
    stop("⚠ Row sums contain NA values. Check input matrix.")
  }

  # Handle samples with zero total count
  zero_rows <- total_counts == 0
  if (any(zero_rows)) {
    warning("⚠ Some samples have zero total read count. CPM will be set to 0 for those samples.")
    total_counts[zero_rows] <- 1
    data[zero_rows, ] <- 0
  }

  # Compute CPM
  cpm_data <- sweep(data, 1, total_counts, FUN = "/") * 1e6

  if (log) {
    cpm_data <- log2(cpm_data + 1)
  }

  return(cpm_data)
}
