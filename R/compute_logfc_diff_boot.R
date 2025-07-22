#' Compute logFC Difference Between Test and Inferred Results (Bootstrapped with Custom Bootstrap Column)
#'
#' Calculates the difference between logFC values from test and inferred datasets
#' for the same features and bootstrap iterations. Allows specification of the
#' bootstrap column and join keys.
#'
#' @param test A data.table with columns `logFC`, a bootstrap column, and join keys.
#' @param infer A data.table with columns `logFC`, a bootstrap column, and join keys.
#' @param join_by A character vector of column names to join by (excluding the bootstrap column).
#' @param boot_col A string giving the name of the bootstrap column (e.g., "bootstrap").
#'
#' @return A data.table with columns: join keys, bootstrap column, `logFC`, and `coef`.
#'
#' @importFrom data.table data.table merge
#' @export
compute_logFC_diff_boot <- function(
  test,
  infer,
  join_by = "feature",
  boot_col = "bootstrap"
) {
  # Validate that the bootstrap column exists in both datasets
  if (!(boot_col %in% names(test)) || !(boot_col %in% names(infer))) {
    stop(sprintf("Both 'test' and 'infer' must contain the bootstrap column '%s'.", boot_col))
  }

  # Validate join keys exist
  missing_test <- setdiff(join_by, names(test))
  missing_infer <- setdiff(join_by, names(infer))
  if (length(missing_test) > 0 || length(missing_infer) > 0) {
    stop("Join columns must exist in both datasets: ",
         paste(c(missing_test, missing_infer), collapse = ", "))
  }

  # Combine full join keys
  full_join_keys <- c(join_by, boot_col)

  # Merge on join keys + bootstrap column
  result <- merge(
    x = test[, c(full_join_keys, list(logFC_test = logFC)), with = FALSE],
    y = infer[, c(full_join_keys, list(logFC_infer = logFC)), with = FALSE],
    by = full_join_keys
  )

  # Compute logFC difference and annotate
  result[, logFC := logFC_test - logFC_infer]
  result[, coef := "test - infer"]

  # Output columns
  output_cols <- c("coef", full_join_keys, "logFC")
  result[, ..output_cols]
}
