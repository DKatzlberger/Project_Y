#' Prepare Expression Data for Selected Features (data.table version)
#'
#' Combines sample metadata with z-scored expression values from a normalized matrix
#' for selected features. Returns long-format data suitable for visualization or analysis.
#'
#' @param meta A data.frame or data.table with sample metadata (samples in rows).
#' @param matrix A numeric matrix of normalized expression values (features x samples).
#' @param features Character vector of feature names to extract (must exist in matrix).
#'
#' @return A data.table in long format with sample metadata, z-scored expression, and feature identifier.
#' @export
prepare_normalized_data <- function(
  meta,
  matrix,
  features
) {
  meta_dt   <- as.data.table(meta)
  meta_dt[, idx := .I]  # create sample index

  exp_list <- lapply(features, function(feat) {
    dt <- copy(meta_dt)
    dt[, zscore := scale(matrix[feat, ])]
    dt[, feature := feat]
    return(dt)
  })

  result <- rbindlist(exp_list)
  return(result[])
}
