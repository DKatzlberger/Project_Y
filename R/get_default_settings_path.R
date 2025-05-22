#' Get and Validate Default Settings File Path
#'
#' Constructs and validates the path to a required YAML settings file.
#' Ensures the file exists and has a valid `.yaml` or `.yml` extension.
#' This function does not load the file — it only returns the path.
#'
#' @param file_name Name of the YAML file (e.g., "default_settings_cross_ancestry_dge.yaml").
#' @param base_dir Directory containing the file, relative to the current working directory.
#'
#' @return A character string with the full, validated file path.
#' @export
get_default_settings_path <- function(
  file_name,
  base_dir
) {
  full_path <- file.path(base_dir, file_name)

  if (!file.exists(full_path)) {
    stop(
      paste0(
        "⚠ Default settings file not found.\n",
        "Expected at: ", normalizePath(full_path, mustWork = FALSE)
      )
    )
  }

  is_yaml_file(full_path)

  return(full_path)
}

