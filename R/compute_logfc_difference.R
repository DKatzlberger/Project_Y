#' Compute logFC Difference Between Test and Inferred Results
#'
#' Calculates the difference between logFC values from test and inferred datasets
#' for the same features. Supports both single observed data and bootstrapped data.
#' If a `bootstrap` column is present, the difference is computed within each iteration.
#'
#' @param test A data.table with columns `feature`, `logFC`, and optionally `bootstrap`.
#' @param infer A data.table with columns `feature`, `logFC`, and optionally `bootstrap`.
#'
#' @return A data.table with columns: `feature`, `logFC`, `coef`, and optionally `bootstrap`.
#'
#' @importFrom data.table data.table merge
#' @export
compute_logfc_difference <- function(
  test,
  infer
) {
  # Determine if bootstrap column is present
  has_bootstrap <- "bootstrap" %in% names(test) && "bootstrap" %in% names(infer)

  # Define merge keys
  merge_keys <- if (has_bootstrap) c("feature", "bootstrap") else "feature"

  # Merge and compute difference
  result <- merge(
    x = test[, .(feature, logFC_test = logFC, bootstrap = if (has_bootstrap) bootstrap else NULL)],
    y = infer[, .(feature, logFC_infer = logFC, bootstrap = if (has_bootstrap) bootstrap else NULL)],
    by = merge_keys
  )

  result[, logFC := logFC_test - logFC_infer]
  result[, coef := "test - infer"]

  # Select final output columns
  output_cols <- c("coef", if (has_bootstrap) "bootstrap", "feature", "logFC")
  result[, ..output_cols]
}
