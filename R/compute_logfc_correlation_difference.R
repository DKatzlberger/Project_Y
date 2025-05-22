#' Compute Correlation Difference Between Train and Target Sets
#'
#' Calculates correlation between logFC values in `train` vs `test` and `train` vs `infer`.
#' Supports both observed data (no `bootstrap`) and bootstrapped data (with `bootstrap`).
#' Returns the difference between correlations in a tidy format.
#'
#' @param train A data.table with columns `feature`, `logFC`, and optionally `bootstrap`.
#' @param test A data.table with columns `feature`, `logFC`, and optionally `bootstrap`.
#' @param infer A data.table with columns `feature`, `logFC`, and optionally `bootstrap`.
#' @param method Correlation method: one of `"pearson"`, `"spearman"`, or `"kendall"`.
#'
#' @return A data.table with columns: `coef`, `method`, `difference`, `cor_train_test`, `cor_train_infer`, and optionally `bootstrap`.
#'
#' @importFrom stats cor
#' @importFrom data.table data.table merge
#' @export
compute_logfc_correlation_difference <- function(
  train,
  test,
  infer,
  method = "pearson"
) {
  has_bootstrap <- "bootstrap" %in% names(test) &&
                   "bootstrap" %in% names(infer)

  merge_keys <- if (has_bootstrap) c("feature", "bootstrap") else "feature"

  merged_test <- merge(
    x = train[, .(feature, bootstrap = if (has_bootstrap) bootstrap else NULL, logFC_train = logFC)],
    y = test[, .(feature, bootstrap = if (has_bootstrap) bootstrap else NULL, logFC_test = logFC)],
    by = merge_keys
  )

  merged_infer <- merge(
    x = train[, .(feature, bootstrap = if (has_bootstrap) bootstrap else NULL, logFC_train = logFC)],
    y = infer[, .(feature, bootstrap = if (has_bootstrap) bootstrap else NULL, logFC_infer = logFC)],
    by = merge_keys
  )

  if (has_bootstrap) {
    # Per-bootstrap correlation
    test_stats <- merged_test[
      , .(cor_train_test = stats::cor(logFC_train, logFC_test, method = method, use = "pairwise.complete.obs")),
      by = bootstrap
    ]

    infer_stats <- merged_infer[
      , .(cor_train_infer = stats::cor(logFC_train, logFC_infer, method = method, use = "pairwise.complete.obs")),
      by = bootstrap
    ]

    result <- merge(test_stats, infer_stats, by = "bootstrap")
    result[, difference := cor_train_test - cor_train_infer]
    result[, `:=`(method = method, coef = "train - test")]

    result[, .(coef, bootstrap, method, cor_train_test, cor_train_infer, difference)]

  } else {
    # Single observed correlation
    cor_test <- stats::cor(
      x      = merged_test$logFC_train,
      y      = merged_test$logFC_test,
      method = method,
      use    = "pairwise.complete.obs"
    )

    cor_infer <- stats::cor(
      x      = merged_infer$logFC_train,
      y      = merged_infer$logFC_infer,
      method = method,
      use    = "pairwise.complete.obs"
    )

    data.table::data.table(
      coef            = "train - test",
      method          = method,
      train_test  = cor_test,
      train_infer = cor_infer,
      difference      = cor_test - cor_infer
    )
  }
}
