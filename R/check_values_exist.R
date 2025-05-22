#' Check That Column Contains Expected Values
#'
#' Validates that all values in a provided vector are present in a data frame column.
#' Stops with an error if any expected values are missing.
#'
#' @param df A data.frame-like object.
#' @param column A character string with the column name to check.
#' @param expected A character vector of expected values.
#'
#' @return Invisible TRUE if all values exist. Otherwise throws an error.
#' @export
check_values_exist <- function(
  df,
  column,
  expected
) {
  if (!column %in% colnames(df)) {
    stop(
      paste0("⚠ Column not found in data: '", column, "'"),
      call. = FALSE
    )
  }

  actual <- unique(df[[column]])
  missing <- setdiff(expected, actual)

  if (length(missing) > 0) {
    stop(
      paste0(
        "⚠ The following expected values were not found in column '", column, "': ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
