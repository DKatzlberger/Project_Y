#' Load and Merge Default and User YAML Settings
#'
#' Loads two YAML files (default and user) and merges them into a single list.
#' User-defined settings override the defaults. Merging is recursive. Adds
#' metadata fields `date` and `id` to the returned settings list.
#'
#' @param default_path File path to the default settings YAML file.
#' @param user_path File path to the user-provided settings YAML file.
#'
#' @return A named list of merged settings with metadata.
#' @export
load_and_merge_settings <- function(
  default_path,
  user_path
) {
  is_yaml_file(
    file_path = default_path
  )

  is_yaml_file(
    file_path = user_path
  )

  default <- yaml::yaml.load_file(
    input = default_path
  )

  user <- yaml::yaml.load_file(
    input = user_path
  )

  setup <- utils::modifyList(
    x         = default,
    val       = user,
    keep.null = TRUE
  )

  setup$date <- format(
    x      = as.POSIXlt(Sys.time(), tz = "GMT"),
    format = "%Y-%m-%d %H:%M:%S"
  )

  setup$id <- toupper(
    substr(
      x     = uuid::UUIDgenerate(),
      start = 1,
      stop  = 10
    )
  )

  return(setup)
}
