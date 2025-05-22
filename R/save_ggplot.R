#' Save a ggplot Object to File with Optional Legend Suppression
#'
#' Safely saves a ggplot object to disk with optional removal of the legend.
#' Includes robust input validation and error handling.
#'
#' @param plot A ggplot object to be saved.
#' @param save_path String. File path where the plot will be saved.
#' @param width Numeric. Width of the plot in inches. Default is 5.
#' @param height Numeric. Height of the plot in inches. Default is 5.
#' @param with_legend Logical. Whether to include the legend. Default is TRUE.
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return Invisibly returns the `save_path` if successful, otherwise `NULL`.
#' @export
save_ggplot <- function(
  plot,
  save_path,
  width       = 5,
  height      = 5,
  with_legend = TRUE,
  ...
) {
  # Validate input: plot
  if (!inherits(plot, "gg")) {
    warning("`plot` must be a valid ggplot object. Plot was not saved.")
    return(NULL)
  }

  # Validate input: save_path
  if (!is.character(save_path) || length(save_path) != 1 || nchar(save_path) == 0) {
    warning("`save_path` must be a non-empty character string. Plot was not saved.")
    return(NULL)
  }

  # Validate input: dimensions
  if (!is.numeric(width) || width <= 0 || !is.numeric(height) || height <= 0) {
    warning("`width` and `height` must be positive numeric values. Plot was not saved.")
    return(NULL)
  }

  # Validate input: legend toggle
  if (!is.logical(with_legend) || length(with_legend) != 1) {
    warning("`with_legend` must be TRUE or FALSE. Plot was not saved.")
    return(NULL)
  }

  # Apply legend removal if needed
  if (!with_legend) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }

  # Try to save the plot
  tryCatch({
    ggplot2::ggsave(
      filename = save_path,
      plot     = plot,
      width    = width,
      height   = height,
      ...
    )
    invisible(save_path)
  }, error = function(e) {
    warning("Failed to save plot to '", save_path, "': ", conditionMessage(e))
    return(NULL)
  })
}
