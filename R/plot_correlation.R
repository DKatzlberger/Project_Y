#' Plot Correlation Bar Plot with Error Bars, Annotations, and Facets
#'
#' @param data A data.frame with correlations and metadata
#' @param x_col Character. Column name for x-axis (e.g. "phase" or "method")
#' @param y_col Character. Column name for y-axis (e.g. "correlation")
#' @param facet_col Character. Column to facet by (e.g. "method")
#' @param pval_col Character. Column name of p-value to annotate (e.g. "p_param")
#' @param ci_lower_col Character. Column name for lower bound of CI
#' @param ci_upper_col Character. Column name for upper bound of CI
#' @param point_size Numeric. Controls text and bar width (default = 0.5)
#' @param title Character. Optional plot title
#'
#' @return A ggplot2 object
#' @export
plot_correlation <- function(
  data,
  x_col,
  y_col,
  facet_col,
  pval_col,
  ci_lower_col = NULL,
  ci_upper_col = NULL,
  point_size = 0.5,
  title = NULL
) {
  library(ggplot2)
  library(rlang)

  facet_formula <- stats::as.formula(paste("~", facet_col))

  # Annotation data (one row per facet level)
  annotation_data <- aggregate(
    data[[pval_col]],
    by = list(data[[facet_col]]),
    FUN = function(p) signif(min(p), 3)
  )
  names(annotation_data) <- c(facet_col, "p_value")

  # Base plot
  p <- ggplot(
    data,
    aes(
      x = !!sym(x_col),
      y = !!sym(y_col)
    )
  ) +
    geom_col(width = point_size) +
    facet_grid(facet_formula)

  # Annotations
  p <- p + geom_text(
    data = annotation_data,
    mapping = aes(
      x = 1,
      y = Inf,
      label = paste0("p = ", p_value)
    ),
    hjust = -0.1,
    vjust = 1.1,
    size = point_size * 3,
    inherit.aes = FALSE
  )

  # Optionally add CI error bars
  if (!is.null(ci_lower_col) && !is.null(ci_upper_col)) {
    p <- p + geom_errorbar(
      aes(
        ymin = !!sym(ci_lower_col),
        ymax = !!sym(ci_upper_col)
      ),
      width = 0.2
    )
  }

  # Final styling
  p <- p +
    labs(
      title = title,
      x = NULL,
      y = y_col
    ) +
    theme_nature_fonts(base_size = point_size * 10) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend(base_size = point_size * 10) +
    theme(legend.position = "none")

  return(p)
}
