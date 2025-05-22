#' Calculate TMM-like Normalization Factors
#'
#' Computes normalization factors using a trimmed mean of M-values approach
#' (similar to edgeR's TMM) for count-based data such as RNA-seq.
#'
#' @param counts         A numeric matrix of raw counts (samples x features).
#' @param logratio_trim  Numeric. Trim proportion for log-ratio values. Default is 0.3.
#' @param abs_expr_trim  Numeric. Trim proportion for absolute expression. Default is 0.05.
#'
#' @return A numeric vector of normalization factors (length = number of samples).
#' @export
calculate_tmm_norm_factors <- function(
  counts,
  logratio_trim = 0.3,
  abs_expr_trim = 0.05
) {

  counts     <- t(counts)  # genes x samples
  lib_sizes  <- colSums(counts)
  ref_index  <- which.min(abs(lib_sizes - median(lib_sizes)))
  ref        <- counts[, ref_index]
  norm_factors <- numeric(ncol(counts))

  for (i in seq_len(ncol(counts))) {
    samp <- counts[, i]

    # Keep genes with non-zero counts in both sample and reference
    keep <- samp > 0 & ref > 0
    x    <- samp[keep]
    y    <- ref[keep]

    if (length(x) < 10) {
      norm_factors[i] <- 1
      next
    }

    M <- log2(x / y)
    A <- 0.5 * log2(x * y)

    m_lower <- quantile(M, logratio_trim / 2)
    m_upper <- quantile(M, 1 - logratio_trim / 2)
    a_lower <- quantile(A, abs_expr_trim / 2)
    a_upper <- quantile(A, 1 - abs_expr_trim / 2)

    keep_trimmed <- (M > m_lower & M < m_upper) & (A > a_lower & A < a_upper)

    if (sum(keep_trimmed) < 10) {
      norm_factors[i] <- 1
    } else {
      norm_factors[i] <- exp(mean(M[keep_trimmed]))
    }
  }

  # Normalize to geometric mean = 1
  gm <- exp(mean(log(norm_factors[norm_factors > 0])))
  norm_factors <- norm_factors / gm

  if (any(norm_factors == 0 | is.na(norm_factors) | is.infinite(norm_factors))) {
    warning("Some normalization factors are zero, NA, or infinite. Check input data.")
  }

  return(norm_factors)
}
