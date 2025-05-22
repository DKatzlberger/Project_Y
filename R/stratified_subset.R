#' Generate Stratified Subset Based on Class Distribution
#'
#' Samples a stratified subset from a data frame using the class distribution
#' provided in `strata`, ensuring reproducibility via a seed.
#'
#' @param data A data frame (e.g., `adata$obs`) containing metadata.
#' @param strata A named integer vector or table with class names and sample counts.
#' @param output_column Character. Column name in `data` representing the class label.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list with `test_idx` and `train_idx`, containing row names.
#' @export
stratified_subset <- function(
  data,
  strata,
  output_column,
  seed
) {
  set.seed(seed)

  test_idx <- character()

  for (stratum in rownames(strata)) {
    n_required    <- strata[stratum]
    stratum_data  <- data[data[[output_column]] == stratum, ]

    if (nrow(stratum_data) >= n_required) {
      sampled_rows <- sample(rownames(stratum_data), n_required)
      test_idx     <- c(test_idx, sampled_rows)
    } else {
      stop(
        paste("Class", stratum, "has fewer rows than required for sampling.")
      )
    }
  }

  train_idx <- rownames(data)[!rownames(data) %in% test_idx]

  return(list(
    test_idx  = test_idx,
    train_idx = train_idx
  ))
}
