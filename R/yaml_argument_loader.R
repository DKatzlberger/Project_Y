#' Load YAML File from Command Line or Use Default (Interactive Mode)
#'
#' Checks whether a YAML file path was passed as a command-line argument.
#' If not, it falls back to a default path and optionally validates it.
#'
#' @param default_yaml A fallback YAML file path for interactive mode.
#' @param example_description A short text label used in interactive mode logging.
#' @param validator A function to validate the file (default: `is_yaml_file()` from utils).
#'
#' @return A string with the path to the selected YAML file.
#' @export
yaml_argument_loader <- function(
  default_yaml        = "example_settings_cross_ancestry_dge.yaml",
  validator           = is_yaml_file
) {
  args <- base::commandArgs(trailingOnly = TRUE)

  if (length(args) > 0) {
    yaml_path <- args[1]
    validator(yaml_path)
  } else {
    base::cat("Running in interactive mode for development.\n")
    base::cat(sprintf("Using example YAML: \"%s\".\n", default_yaml))
    yaml_path <- default_yaml
    validator(yaml_path)
  }

  return(yaml_path)
}