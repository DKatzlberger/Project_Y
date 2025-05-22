#' Stratify Dataset by Ancestry and Class
#'
#' Splits a dataset into training, testing, and inference subsets based on ancestry,
#' using stratified sampling by class. Optionally saves a stratification QC plot.
#'
#' @param adata AnnData-like object with `$obs` for metadata and `$X` for data.
#' @param ancestry_column Character. Column name for ancestry labels in `adata$obs`.
#' @param output_column Character. Column name for class labels to stratify on.
#' @param train_ancestry Character vector. Ancestry values used for train/test.
#' @param infer_ancestry Character vector. Ancestry values held out for inference.
#' @param seed Integer. Random seed for reproducibility.
#' @param plot_path Optional. File path to save a stratification QC plot.
#'
#' @return A list with `train_adata`, `test_adata`, and `infer_adata`.
#' @export
stratify_ancestry_subsets <- function(
  adata,
  ancestry_column,
  output_column,
  train_ancestry,
  infer_ancestry,
  seed,
  plot_path = NULL
) {
  # Subset by ancestry
  train_adata_ <- adata[adata$obs[[ancestry_column]] %in% train_ancestry, ]
  infer_adata  <- adata[adata$obs[[ancestry_column]] %in% infer_ancestry, ]

  # Get class distribution from infer ancestry
  strata <- table(infer_adata$obs[[output_column]])

  # Stratify train ancestry into train/test
  indices <- stratified_subset(
    data          = train_adata_$obs,
    strata        = strata,
    output_column = output_column,
    seed          = seed
  )

  stopifnot(length(intersect(indices$train_idx, indices$test_idx)) == 0)

  train_adata <- train_adata_[indices$train_idx, ]
  test_adata  <- train_adata_[indices$test_idx, ]

  # Optional QC plot
  if (!is.null(plot_path)) {
    meta_combined <- dplyr::bind_rows(
      dplyr::mutate(train_adata$obs, Set = "Train"),
      dplyr::mutate(test_adata$obs,  Set = "Test"),
      dplyr::mutate(infer_adata$obs, Set = "Infer")
    )

    p_count <- plot_output_column_count(
      data = meta_combined,
      x    = "Set",
      fill = output_column
    )

    p_prop <- plot_output_column_proportion(
      data = meta_combined,
      x    = "Set",
      fill = output_column
    )

    p_combined <- p_count + p_prop + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")

    # Estimate optimal plot size
    size <- estimate_plot_size(p_combined)

    save_ggplot(
      plot      = p_combined,
      save_path = plot_path,
      width     = size$width,
      height    = size$height
    )
  }

  return(list(
    train_adata = train_adata,
    test_adata  = test_adata,
    infer_adata = infer_adata
  ))
}
