#' Plot Null Distribution with Observed Statistic, Fitted Normal Curve, and P-values
#'
#' @param observed_dt A data.table or data.frame containing observed statistics.
#' @param bootstrap_dt A data.table or data.frame containing bootstrap samples (null distribution).
#' @param value_col Character. Column name of the numeric test statistic (e.g., "logFC").
#' @param by_col Optional character. Column name for faceting/grouping (e.g., "feature"). Use NULL for global.
#' @param features Optional character vector. A subset of `by_col` values to include in the plot.
#' @param point_size Numeric. Controls the size of lines and text scaling.
#'
#' @return A ggplot object.
#' @export
plot_null_distribution <- function(
  observed_dt,
  bootstrap_dt,
  value_col  = "logFC",
  by_col     = NULL,
  features   = NULL,
  point_size = 0.5
) {
  data.table::setDT(observed_dt)
  data.table::setDT(bootstrap_dt)

  if (!is.null(by_col)) {
    observed_dt   <- observed_dt[get(by_col) %in% features]
    bootstrap_dt  <- bootstrap_dt[get(by_col) %in% features]

    obs_vals  <- observed_dt[, .(observed = get(value_col)), by = by_col]
    boot_vals <- bootstrap_dt[, .(value = get(value_col)), by = by_col]
    data      <- merge(boot_vals, obs_vals, by = by_col)

    # Compute p-values
    data[, empirical_p := (1 + sum(abs(value) >= abs(observed))) / (.N + 1), by = by_col]
    data[, mu := mean(value), by = by_col]
    data[, sigma := sd(value), by = by_col]
    data[, z := abs((observed - mu) / sigma)]
    data[, normal_p := 2 * stats::pnorm(-z)]

    # Create annotation data
    label_data <- data[, .SD[1], by = by_col]
    label_data[, label := sprintf("p_emp = %.3g\np_norm = %.3g", empirical_p, normal_p)]

    plot <- ggplot2::ggplot(
      data    = data,
      mapping = ggplot2::aes(x = value)
    ) +
      ggplot2::geom_histogram(
        mapping = ggplot2::aes(y = ggplot2::after_stat(density)),
        bins    = 30
      ) +
      ggplot2::geom_vline(
        mapping   = ggplot2::aes(xintercept = observed),
        color     = "orange",
        linewidth = point_size / 1.5
      ) +
      ggplot2::stat_function(
        fun       = stats::dnorm,
        args      = list(mean = mean(data$value), sd = sd(data$value)),
        color     = "blue",
        linewidth = point_size / 1.5
      ) +
      ggplot2::geom_text(
        data    = label_data,
        mapping = ggplot2::aes(
          x     = -Inf,
          y     = Inf,
          label = label
        ),
        hjust = -0.05,
        vjust = 1.2,
        size  = point_size * 3,
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", by_col))
      ) +
      ggplot2::labs(
        x = value_col,
        y = "Density"
      ) +
      theme_nature_fonts(
        base_size = point_size * 10
      ) +
      theme_white_background() +
      theme_white_strip() +
      theme_small_legend(
        base_size = point_size * 10
      )

    return(plot)
  }

  # Global (non-faceted) version
  t_obs  <- observed_dt[[value_col]]
  t_null <- bootstrap_dt[[value_col]]
  mu     <- mean(t_null)
  sigma  <- sd(t_null)

  empirical_p <- (1 + sum(abs(t_null) >= abs(t_obs))) / (length(t_null) + 1)
  z           <- abs((t_obs - mu) / sigma)
  normal_p    <- 2 * stats::pnorm(-z)
  label_text  <- sprintf("p_emp = %.3g\np_norm = %.3g", empirical_p, normal_p)

  data <- data.table::data.table(value = t_null)

  plot <- ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(x = value)
  ) +
    ggplot2::geom_histogram(
      mapping = ggplot2::aes(y = ggplot2::after_stat(density)),
      bins    = 30
    ) +
    ggplot2::geom_vline(
      xintercept = t_obs,
      color      = "orange",
      linewidth  = point_size / 1.5
    ) +
    ggplot2::stat_function(
      fun       = stats::dnorm,
      args      = list(mean = mu, sd = sigma),
      color     = "blue",
      linewidth = point_size / 1.5
    ) +
    ggplot2::geom_text(
      data = data.frame(label = label_text),
      mapping = ggplot2::aes(
        x     = -Inf,
        y     = Inf,
        label = label
      ),
      hjust = -0.05,
      vjust = 1.2,
      size  = point_size * 3,
      inherit.aes = FALSE
    ) + 
    ggplot2::labs(
      x = value_col,
      y = "Density"
    ) +
    theme_nature_fonts(
      base_size = point_size * 10
    ) +
    theme_white_background() +
    theme_white_strip() +
    theme_small_legend(
      base_size = point_size * 10
    )

  return(plot)
}
