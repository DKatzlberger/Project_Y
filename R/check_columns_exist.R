#' Check if Required Columns Exist in DataFrame
#'
#' @param df A data.frame-like object (e.g. adata$obs)
#' @param required_cols A character vector of required column names
#'
#' @return TRUE (invisible); stops with error if any are missing.
#' @export
check_columns_exist <- function(
  df,
  required_cols
) {
  missing <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing) > 0) {
    stop(
      paste0(
        "⚠ Missing required columns in data: ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}