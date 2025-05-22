#' Calculate Log Fold Change Between Two Groups
#'
#' Computes the log fold change (logFC) between two levels of a grouping factor
#' across all features in a matrix, optionally applying normalization (e.g., voom).
#'
#' @param matrix A numeric matrix of expression values with samples in rows and features (e.g., genes) in columns.
#' @param meta A data frame of sample metadata. Row names must match those of `matrix`.
#' @param group_column A character string specifying the name of the grouping column in `meta`.
#'   Must be a factor with exactly two levels.
#' @param normalization An optional normalization function. Must accept a transposed matrix
#'   (features x samples) and a design matrix, and return either a matrix or
#'   a list with an `E` element containing normalized values.
#'
#' @return A `data.table` with one row per feature, containing the columns:
#'   \describe{
#'     \item{coef}{A string representing the group contrast (e.g., "A - B")}
#'     \item{feature}{The feature (e.g., gene) name}
#'     \item{logFC}{The estimated log fold change between groups}
#'   }
#'
#' @importFrom matrixStats colMeans2
#' @importFrom data.table data.table
#' @export
calculate_logfc <- function(
  matrix,
  meta,
  group_column,
  normalization = NULL
) {
  stopifnot(
    all(rownames(meta) == rownames(matrix))
  )

  groups <- meta[[group_column]]

  if (!is.factor(groups)) {
    stop("`group_column` must be a factor.")
  }

  if (length(levels(groups)) != 2) {
    stop("`group_column` must have exactly 2 levels.")
  }

  group_1 <- levels(groups)[1]
  group_2 <- levels(groups)[2]

  group_1_idx <- groups == group_1
  group_2_idx <- groups == group_2

  if (!is.null(normalization)) {
    design_formula <- stats::as.formula(paste("~0 +", group_column))
    design_matrix <- stats::model.matrix(design_formula, data = meta)

    rownames(design_matrix) <- rownames(meta)

    normalized <- normalization(
      t(matrix),
      design_matrix
    )

    expression_matrix <- if (is.list(normalized) && !is.null(normalized$E)) {
      t(normalized$E)
    } else {
      t(normalized)
    }
  } else {
    expression_matrix <- matrix
  }

  mean_1 <- matrixStats::colMeans2(
    expression_matrix[group_1_idx, , drop = FALSE]
  )

  mean_2 <- matrixStats::colMeans2(
    expression_matrix[group_2_idx, , drop = FALSE]
  )

  logfc <- mean_1 - mean_2

  data.table::data.table(
    coef    = paste(group_1, "-", group_2),
    feature = colnames(matrix),
    logFC   = logfc
  )
}
