#' Validate That a File is an H5AD File
#'
#' Checks that the provided file path exists and has a `.h5ad` extension.
#' Throws an error if the file is missing or not in the correct format.
#'
#' @param file_path Path to the file to check.
#'
#' @return Invisibly returns TRUE if valid; otherwise throws an error.
#' @export
is_h5ad_file <- function(
  file_path
) {
  if (!file.exists(file_path)) {
    stop(
      paste0("⚠ File does not exist: ", file_path),
      call. = FALSE
    )
  }

  ext <- tolower(tools::file_ext(file_path))

  if (ext != "h5ad") {
    stop(
      paste0("⚠ File is not an .h5ad file: ", file_path),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
