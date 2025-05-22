#' Filter Features by Signal and Sample Coverage
#'
#' Filters features (e.g., genes) based on total signal and the number of samples
#' in which they are expressed. Useful for excluding low-expression features.
#'
#' @param data A numeric matrix of expression values (samples x features).
#' @param min_signal Numeric. Minimum total signal across all samples. Default is 10.
#' @param max_signal Optional numeric. Maximum total signal (ignored if NULL).
#' @param min_samples_ratio Numeric. Minimum fraction of samples in which a feature
#'        must be expressed (> 0). Default is 0.5 (i.e., 50% of samples).
#'
#' @return A logical vector indicating which features to keep.
#' @export
filter_by_signal <- function(
  data,
  min_signal        = 10,
  max_signal        = NULL,
  min_samples_ratio = 0.5
) {
  if (!is.matrix(data)) {
    stop("⚠ `data` must be a numeric matrix.")
  }

  if (!is.numeric(min_signal) || length(min_signal) != 1) {
    stop("⚠ `min_signal` must be a single numeric value.")
  }

  if (!is.null(max_signal) && (!is.numeric(max_signal) || length(max_signal) != 1)) {
    stop("⚠ `max_signal` must be NULL or a single numeric value.")
  }

  if (!is.numeric(min_samples_ratio) || min_samples_ratio <= 0 || min_samples_ratio > 1) {
    stop("⚠ `min_samples_ratio` must be a numeric value between 0 and 1.")
  }

  num_samples   <- nrow(data)
  min_samples   <- floor(min_samples_ratio * num_samples)
  gene_signal   <- colSums(data, na.rm = TRUE)
  gene_detected <- colSums(data > 0, na.rm = TRUE)

  if (is.null(max_signal)) {
    keep <- (gene_signal >= min_signal) &
            (gene_detected >= min_samples)
  } else {
    keep <- (gene_signal >= min_signal) &
            (gene_signal <= max_signal) &
            (gene_detected >= min_samples)
  }

  return(keep)
}
