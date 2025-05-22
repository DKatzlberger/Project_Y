#' Estimate Optimal Plot Size from a ggplot Object
#'
#' Infers an appropriate width and height (in inches) for a ggplot object,
#' based on facets, axis labels, legend, and base sizing heuristics.
#'
#' @param plot A ggplot object.
#' @param base_unit Numeric. Base size (in inches) for a single facet panel. Default is 2.
#' @param text_scaling Numeric. Multiplier to adjust for axis label and legend size. Default is 1.
#' @param legend_padding Numeric. Extra inches to add if a legend is present. Default is 1.
#'
#' @return A named list with elements `width` and `height`.
#' @export
estimate_plot_size <- function(
  plot,
  base_unit      = 2,
  text_scaling   = 1,
  legend_padding = 1
) {
  if (!inherits(plot, "gg")) {
    stop("`plot` must be a ggplot object.")
  }

  # Build the ggplot object to access layout
  built <- ggplot2::ggplot_build(plot)
  layout <- built$layout

  # Get facet dimensions
  nrow_facets <- layout$facet$params$nrow %||% 1
  ncol_facets <- layout$facet$params$ncol %||% 1

  # Estimate axis label length scaling (based on string width)
  axis_text_x <- built$plot$labels$x %||% ""
  axis_text_y <- built$plot$labels$y %||% ""
  x_len <- nchar(axis_text_x)
  y_len <- nchar(axis_text_y)
  axis_scaling <- 1 + 0.05 * max(x_len, y_len)

  # Estimate legend impact
  has_legend <- !inherits(built$plot$theme$legend.position, "element_blank") &&
    (!is.null(built$plot$theme$legend.position) && built$plot$theme$legend.position != "none")

  # Count discrete legend items if fill or color aesthetics used
  n_legend_items <- 0
  mapped_aes <- names(built$plot$mapping)
  discrete_vars <- c("fill", "color", "colour")

  for (aes in discrete_vars) {
    if (aes %in% mapped_aes) {
      var_name <- rlang::as_name(built$plot$mapping[[aes]])
      data_col <- built$data[[1]][[var_name]]
      if (is.factor(data_col) || is.character(data_col)) {
        n_legend_items <- max(n_legend_items, length(unique(data_col)))
      }
    }
  }

  legend_adjustment <- if (has_legend && n_legend_items > 0) {
    legend_padding + (n_legend_items / 10)
  } else {
    0
  }

  # Final estimated size
  width <- ncol_facets * base_unit * axis_scaling * text_scaling + legend_adjustment
  height <- nrow_facets * base_unit * axis_scaling * text_scaling

  list(
    width  = round(width, 2),
    height = round(height, 2)
  )
}
