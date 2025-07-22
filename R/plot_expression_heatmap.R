#' Plot and Save Gene Expression Heatmap with Metadata Annotations
#'
#' Plots a heatmap of gene expression z-scores with annotations for ancestry and output subtype.
#'
#' @param zscores A data frame or data.table with at least 'idx', 'feature', 'zscore', and metadata columns.
#' @param ancestry_column Name of the column to use for ancestry annotation.
#' @param output_column Name of the column to use for output subtype annotation.
#' @param file_path Optional file path to save the plot (e.g., "heatmap.pdf", "heatmap.png").
#' @param cluster_rows, cluster_cols Logical; whether to cluster rows or columns. Default: TRUE.
#' @param show_row_names Logical; show row (gene) names. Default: TRUE.
#' @param color_scale A color function for heatmap values. Default is blue-white-red for z-scores.
#' @param width, height Width and height in inches (if saving). Optional.
#' @param units Units for plot size when saving. Default is "in".
#' @param res Resolution in dpi for PNG/TIFF output. Default is 300.
#' @param point_size Numeric. Font size scaling. Default is 0.5 (~5pt).
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap`.
#'
#' @return Invisibly returns the z-score matrix used for the heatmap.
#' @export
plot_expression_heatmap <- function(
  zscores,
  ancestry_column,
  output_column,
  file_path      = NULL,
  color_scale    = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  width          = NULL,
  height         = NULL,
  units          = "in",
  res            = 300,
  point_size     = 0.5,
  ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  zscores <- data.table::as.data.table(zscores)

  # Validate required columns
  required_cols <- c("idx", "feature", "zscore", ancestry_column, output_column)
  missing_cols <- setdiff(required_cols, names(zscores))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in zscores: ", paste(missing_cols, collapse = ", "))
  }

  # Wide format: features x samples
  expr_wide             <- data.table::dcast(zscores, feature ~ idx, value.var = "zscore")
  expr_matrix           <- as.matrix(expr_wide[, -1, with = FALSE])
  rownames(expr_matrix) <- expr_wide$feature

  # Metadata
  annot_df <- unique(zscores[, .(idx, ancestry = get(ancestry_column), output = get(output_column))])
  annot_df <- annot_df[match(colnames(expr_matrix), annot_df$idx)]

  if (anyNA(annot_df$idx)) {
    stop("Some sample IDs in the expression matrix were not matched with metadata.")
  }

  # Sort samples by ancestry, then output
  annot_df <- annot_df[order(annot_df$ancestry, annot_df$output)]
  expr_matrix <- expr_matrix[, annot_df$idx]

  # Annotation dataframe
  annotation_df <- setNames(
    data.frame(
      annot_df$ancestry,
      annot_df$output,
      row.names = annot_df$idx,
      check.names = FALSE
    ),
    c(ancestry_column, output_column)
  )

  # Color mapping
  ancestry_levels <- unique(annotation_df[[ancestry_column]])
  output_levels   <- unique(annotation_df[[output_column]])

  ancestry_colors <- setNames(c("navy", "goldenrod"), ancestry_levels)
  output_colors   <- setNames(c("purple", "#027c58"), output_levels)

  # Annotations
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = setNames(
      list(ancestry_colors, output_colors),
      c(ancestry_column, output_column)
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(direction = "horizontal")
  )

  # Define column split (factor with custom levels)
  col_split <- factor(
    interaction(annotation_df[[ancestry_column]], annotation_df[[output_column]], drop = TRUE),
    levels = unique(interaction(annotation_df[[ancestry_column]], annotation_df[[output_column]]))
  )

  # Build heatmap
  heatmap_obj <- ComplexHeatmap::Heatmap(
    expr_matrix,
    name              = "Z-score",
    col               = color_scale,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    row_names_side    = "left",
    column_split      = col_split,
    column_title      = NULL,  
    top_annotation    = top_anno,
    heatmap_legend_param = list(direction = "horizontal"),
    ...
  )

  # Optional: Save to file
  if (!is.null(file_path)) {
    ext <- tools::file_ext(file_path)
    switch(
      tolower(ext),
      pdf  = grDevices::pdf(file_path, width = width %||% 7, height = height %||% 7),
      png  = grDevices::png(file_path, width = width %||% 7, height = height %||% 7, units = units, res = res),
      tiff = grDevices::tiff(file_path, width = width %||% 7, height = height %||% 7, units = units, res = res),
      stop("Unsupported file extension: ", ext)
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Draw heatmap
  ComplexHeatmap::draw(
    heatmap_obj,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "bottom",
    show_heatmap_legend    = TRUE,
    show_annotation_legend = TRUE,
    merge_legend           = TRUE
  )

  invisible(expr_matrix)
}
