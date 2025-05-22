#' Calculate empirical and parametric p-values from bootstrap distribution
#'
#' @description Computes two-sided empirical and parametric (normal-theory) p-values
#' comparing observed statistics to a null distribution obtained via bootstrapping.
#' Optionally creates and saves a plot of the null distribution(s).
#'
#' @param observed_dt A data.table containing observed test statistics.
#' @param bootstrap_dt A data.table containing bootstrap samples (null distribution).
#' @param value_col String. Name of the numeric column to compare (e.g., "difference" or "logFC").
#' @param by_cols String or character vector. Grouping columns (e.g., "feature"). Ignored if `is_global = TRUE`.
#' @param is_global Logical. If TRUE, treats the data as a single global test (no grouping).
#' @param plot_path Optional string. If provided, a plot of the null distribution(s) is saved to this file path.
#' @param features Optional character vector. A subset of features to include in the plot (only applies when grouped).
#'
#' @return A data.table with columns: grouping variables (if any), `method`, observed column, `cor_train_test`, `cor_train_infer`, `p_param`, and `p_emp`.
#' @export
calculate_pvalues <- function(
  observed_dt,
  bootstrap_dt,
  value_col  = "logFC",
  by_cols    = "feature",
  is_global  = FALSE,
  plot_path  = NULL,
  features   = NULL
) {
  data.table::setDT(observed_dt)
  data.table::setDT(bootstrap_dt)
  observed_dt  <- data.table::copy(observed_dt)
  bootstrap_dt <- data.table::copy(bootstrap_dt)

  if (is_global) {
    t_obs  <- observed_dt[[value_col]]
    t_null <- bootstrap_dt[[value_col]]
    b      <- length(t_null)

    p_emp   <- (1 + sum(abs(t_null) >= abs(t_obs))) / (1 + b)
    mu      <- mean(t_null)
    sigma   <- sd(t_null)
    z       <- abs((t_obs - mu) / sigma)
    p_param <- 2 * stats::pnorm(-z)

    result <- data.table::data.table(
      method          = observed_dt$method,
      observed_value  = t_obs,
      train_test      = observed_dt$train_test,
      train_infer     = observed_dt$train_infer,
      p_param         = p_param,
      p_emp           = p_emp
    )

    setnames(result, "observed_value", value_col)

    if (!is.null(plot_path)) {
      p <- plot_null_distribution(
        observed_dt  = observed_dt,
        bootstrap_dt = bootstrap_dt,
        value_col    = value_col,
        by_col       = NULL
      )

      size <- estimate_plot_size(p)

      save_ggplot(
        save_path = plot_path,
        plot      = p,
        width     = size$width,
        height    = size$height
      )
    }

    return(result[])
  }

  # Grouped case
  observed_dt[, observed := .SD[[value_col]], .SDcols = value_col]
  bootstrap_dt[, value := .SD[[value_col]], .SDcols = value_col]

  merged_dt <- merge(
    bootstrap_dt[, c(by_cols, "value"), with = FALSE],
    observed_dt[, c(by_cols, "observed"), with = FALSE],
    by = by_cols
  )

  result <- merged_dt[, .(
    observed = unique(observed),
    p_emp    = (1 + sum(abs(value) >= abs(observed))) / (1 + .N),
    p_param  = {
      mu    <- mean(value)
      sigma <- sd(value)
      z     <- abs((unique(observed) - mu) / sigma)
      2 * stats::pnorm(-z)
    }
  ), by = by_cols]

  setnames(result, "observed", value_col)

  if (!is.null(plot_path)) {
    plot_features <- if (is.null(features)) result[[by_cols]] else features

    p <- plot_null_distribution(
      observed_dt  = observed_dt,
      bootstrap_dt = bootstrap_dt,
      value_col    = value_col,
      by_col       = by_cols,
      features     = plot_features
    )

    size <- estimate_plot_size(p)

    save_ggplot(
      save_path = plot_path,
      plot      = p,
      width     = size$width,
      height    = size$height
    )
  }

  return(result[])
}
