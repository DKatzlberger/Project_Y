#' Check for Required Settings Keys
#'
#' Ensures that all required keys are present in a settings list.
#' If any keys are missing, the function will stop with an informative error.
#'
#' @param settings A named list of loaded and merged settings.
#' @param required_keys A character vector of required top-level keys.
#'
#' @return Invisible TRUE if validation passes; otherwise throws an error.
#' @export
check_required_settings <- function(
  settings,
  required_keys
) {
  missing_keys <- setdiff(
    required_keys,
    names(settings)
  )

  if (length(missing_keys) > 0) {
    msg <- paste0(
      "⚠ Missing required settings:\n",
      paste("  -", missing_keys, collapse = "\n"),
      "\n\nPlease add them to your YAML configuration."
    )
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}