#' Normalize Count Matrix Using Voom (No edgeR)
#'
#' Applies voom normalization from the limma package to a count matrix.
#'
#' @param count_matrix A numeric matrix of raw counts (genes x samples).
#' @param design_matrix A model matrix for the experimental design (samples x covariates).
#'
#' @return An `EList` object from `limma::voom()`, including normalized logCPM values in `E`.
#' @export
voom_normalization <- function(
  count_matrix,
  design_matrix
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("The 'limma' package must be installed to use voom normalization.")
  }

  # Compute library sizes and normalization factors manually
  lib_sizes <- colSums(count_matrix)
  norm_factors <- lib_sizes / mean(lib_sizes)

  # Apply voom using manually calculated normalization factors
  voom_res <- limma::voom(counts   = count_matrix,
                          design   = design_matrix,
                          lib.size = lib_sizes * norm_factors,
                          plot     = FALSE
                          )

  return(voom_res)
}