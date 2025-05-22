#' Extract differential expression results from limma model
#'
#' @param limma_fit A limma MArrayLM fit object (from `lmFit()` or `contrasts.fit()`).
#' @param feature_column Name of the column that stores rownames (gene IDs). Default is "feature".
#' @param exclude_intercept Logical; if TRUE (default), removes rows with `(Intercept)` coefficient.
#'
#' @return A `data.table` with all topTable results, stacked and labeled by coefficient. `coef` is the first column.
#' @export
#'
#' @examples
#' extract_top_table(limma_fit, feature_column = "gene", exclude_intercept = TRUE)
extract_top_table <- function(
  limma_fit,
  feature_column,
  exclude_intercept = TRUE
) {
  stopifnot(inherits(limma_fit, "MArrayLM"))

  results_list <- lapply(
    colnames(coef(limma_fit)),
    function(coef_name) {
      tbl <- limma::topTable(
        fit = limma_fit,
        coef = coef_name,
        number = Inf
      )

      tbl_dt <- data.table::as.data.table(tbl, keep.rownames = feature_column)
      tbl_dt[, coef := coef_name]
      tbl_dt
    }
  )

  results_combined <- data.table::rbindlist(results_list)

  if (exclude_intercept) {
    results_combined <- results_combined[coef != "(Intercept)"]
  }

  # Reorder columns: coef first, then feature_column, then everything else
  setcolorder(
    results_combined,
    c("coef", feature_column, setdiff(colnames(results_combined), c("coef", feature_column)))
  )

  return(results_combined[])
}

