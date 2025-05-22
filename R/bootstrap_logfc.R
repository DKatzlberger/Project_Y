#' Bootstrap Log Fold Change Between Two Groups
#'
#' Computes log fold change estimates via bootstrap sampling between two groups.
#'
#' @param matrix A numeric matrix with samples in rows and features in columns.
#' @param meta A data frame containing sample annotations (must have same rownames as `matrix`).
#' @param group_column The name of the column in `meta` defining two groups (must be a factor).
#' @param size Number of samples to draw in each bootstrap iteration.
#' @param normalization Optional normalization function. Should accept `(t(matrix), design_matrix)`.
#' @param n_iterations Number of bootstrap iterations to run. Default is 1000.
#' @param n_cpus Number of CPU cores to use for parallel execution. Defaults to `floor(n_iterations / 100)`.
#'
#' @return A data.table with logFC values per bootstrap iteration.
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom matrixStats colMeans2
#' @importFrom furrr future_map
#' @importFrom future plan multisession
#' @export
bootstrap_logfc <- function(
  matrix,
  meta,
  group_column,
  size,
  normalization = NULL,
  n_iterations  = 1000,
  n_cpus        = 1
) {
  stopifnot(all(rownames(meta) == rownames(matrix)))

  max_safe_cores <- parallel::detectCores() - 1

  if (is.null(n_cpus)) {
    n_cpus <- max(1, floor(n_iterations / 100))
  }

  if (n_cpus > max_safe_cores) {
    message(
      sprintf("⚠️ Requested %d cores; limiting to %d due to system constraints.", n_cpus, max_safe_cores)
    )
    n_cpus <- max_safe_cores
  }

  message(sprintf("Running %d bootstrap iterations using %d cores.", n_iterations, n_cpus))

  future::plan(future::multisession, workers = n_cpus)

  start_time <- Sys.time()

  results_list <- furrr::future_map(
    .x       = seq_len(n_iterations),
    .options = furrr::furrr_options(seed = TRUE),
    .f       = function(i) {
      tryCatch({
        idx <- sample(seq_len(nrow(matrix)), size = size, replace = TRUE)

        boot_matrix <- matrix[idx, , drop = FALSE]
        boot_meta   <- meta[idx, , drop = FALSE]

        groups <- boot_meta[[group_column]]

        if (!is.factor(groups)) {
          stop("Group column must be a factor.")
        }

        if (length(levels(groups)) != 2) {
          stop("Group factor must have exactly 2 levels.")
        }

        g1 <- levels(groups)[1]
        g2 <- levels(groups)[2]

        if (!is.null(normalization)) {
          rownames(boot_matrix) <- paste0("Sample_", i, "_", seq_len(nrow(boot_matrix)))
          rownames(boot_meta)   <- rownames(boot_matrix)

          design <- stats::model.matrix(
            stats::as.formula(paste("~0 +", group_column)),
            data = boot_meta
          )

          norm <- normalization(t(boot_matrix), design)

          boot_matrix <- if (is.list(norm) && !is.null(norm$E)) {
            t(norm$E)
          } else {
            t(norm)
          }
        }

        g1_idx <- groups == g1
        g2_idx <- groups == g2

        mean_1 <- matrixStats::colMeans2(boot_matrix[g1_idx, , drop = FALSE])
        mean_2 <- matrixStats::colMeans2(boot_matrix[g2_idx, , drop = FALSE])
        logfc  <- mean_1 - mean_2

        data.table::data.table(
          coef      = paste(g1, "-", g2),
          feature   = colnames(matrix),
          logFC     = logfc,
          bootstrap = i
        )
      }, error = function(e) {
        message(sprintf("⚠️ Iteration %d failed: %s", i, e$message))
        NULL
      })
    }
  )

  duration <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
  message(sprintf("Bootstrap completed in %s seconds.", duration))

  data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
}
