#' Validate a YAML File Path
#'
#' Checks whether the provided file path exists and has a valid `.yaml` or `.yml` extension.
#'
#' @param file_path A character string giving the path to the file.
#'
#' @return Returns `TRUE` invisibly if the file is valid; otherwise, throws an error.
#' @export
is_yaml_file <- function(file_path) {
  # Check that the file exists
  if (!base::file.exists(file_path)) {
    stop(sprintf("⚠ The file '%s' does not exist.", file_path))
  }

  # Check the file extension
  extension <- base::tolower(tools::file_ext(file_path))
  if (!(extension %in% c("yml", "yaml"))) {
    stop("⚠ Invalid file extension. Expected a '.yaml' or '.yml' file.")
  }

  # Return TRUE invisibly for piping or conditional use
  return(invisible(TRUE))
}
