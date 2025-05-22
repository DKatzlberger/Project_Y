#' Save Final Settings to YAML
#'
#' Creates the output directory if it doesn't exist and writes
#' the settings list to a YAML file at the given path.
#'
#' @param settings A named list of configuration settings.
#' @param file_path Full path to the YAML file to write.
#'
#' @return The path to the saved file (invisible).
#' @export
save_settings <- function(
  settings,
  file_path
) {
  output_dir <- dirname(file_path)

  if (!dir.exists(output_dir)) {
    dir.create(
      path      = output_dir,
      recursive = TRUE
    )
  }

  yaml::write_yaml(
    x    = settings,
    file = file_path
  )

  invisible(file_path)
}