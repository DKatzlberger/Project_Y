#' Compute logFC Difference Between Test and Inferred Results (Non-Bootstrapped)
#'
#' Calculates the difference between logFC values from test and inferred datasets
#' for the same features. This version assumes no `bootstrap` column is present.
#' Allows custom join keys via the `join_by` argument.
#'
#' @param test A data.table with columns including `logFC` and join keys.
#' @param infer A data.table with columns including `logFC` and join keys.
#' @param join_by A character vector specifying the column(s) to join on (default: "feature").
#'
#' @return A data.table with columns: join keys, `logFC`, and `coef`.
#'
#' @importFrom data.table data.table merge
#' @export
compute_logFC_diff_obs <- function(
  test,
  infer,
  join_by = "feature"
) {
  # Validate no 'bootstrap' column
  if ("bootstrap" %in% names(test) || "bootstrap" %in% names(infer)) {
    stop("This function does not support data with a 'bootstrap' column.")
  }

  # Ensure all join_by columns exist in both inputs
  missing_test <- setdiff(join_by, names(test))
  missing_infer <- setdiff(join_by, names(infer))
  if (length(missing_test) > 0 || length(missing_infer) > 0) {
    stop("Join columns must exist in both datasets: ",
         paste(c(missing_test, missing_infer), collapse = ", "))
  }

  # Merge on join_by
  result <- merge(
    x = test[, c(join_by, list(logFC_test = logFC)), with = FALSE],
    y = infer[, c(join_by, list(logFC_infer = logFC)), with = FALSE],
    by = join_by
  )

  # Compute the difference and label
  result[, logFC := logFC_test - logFC_infer]
  result[, coef := "test - infer"]

  # Output: join keys, coef, and difference
  output_cols <- c("coef", join_by, "logFC")
  result[, ..output_cols]
}
