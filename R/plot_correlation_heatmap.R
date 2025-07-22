#' Plot and Save Clustered Correlation Heatmap of Samples with Metadata Annotations
#'
#' Computes a correlation matrix and plots a heatmap with annotations for ancestry and output.
#'
#' @param expr_data Matrix or data frame: genes x samples (or samples x genes — auto-detected).
#' @param method Correlation method: "pearson", "spearman", or "kendall". Default: "pearson".
#' @param file_path Optional file path to save the plot (e.g., "plot.pdf", "plot.png").
#' @param metadata Optional data frame with sample metadata (rownames must match colnames of expr_data).
#' @param ancestry_column Metadata column name to annotate on both rows and columns (2-level factor).
#' @param output_column Metadata column name to annotate on both rows and columns (2-level factor).
#' @param cluster_rows, cluster_cols Logical; cluster heatmap rows/columns. Default is TRUE.
#' @param show_names Logical; show row/column names. Default is TRUE.
#' @param color_scale A color function for the heatmap body. Default is blue-white-red.
#' @param width, height Width and height in inches (if saving). Optional.
#' @param units Units for plot size if saving. Default is "in".
#' @param res Resolution in dpi for PNG/TIFF output. Default is 300.
#' @param point_size Numeric. Font size scaling for sample names and legends. Default is 0.5 (~5pt).
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap`.
#'
#' @return Invisibly returns the correlation matrix.
#' @export
plot_correlation_heatmap <- function(
  expr_data,
  method           = "pearson",
  file_path        = NULL,
  metadata         = NULL,
  ancestry_column  = NULL,
  output_column    = NULL,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  show_names       = TRUE,
  color_scale      = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  width            = NULL,
  height           = NULL,
  units            = "in",
  res              = 300,
  point_size       = 0.5,
  ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!is.matrix(expr_data) && !is.data.frame(expr_data)) {
    stop("`expr_data` must be a matrix or data frame.")
  }

  expr_data <- as.matrix(expr_data)

  if (nrow(expr_data) < ncol(expr_data)) {
    warning("Input appears to be samples x genes. Transposing to compute sample correlations.")
    expr_data <- t(expr_data)
  }

  cor_mat <- cor(expr_data, method = method)
  sample_ids <- colnames(expr_data)

  top_anno <- NULL
  left_anno <- NULL

  if (!is.null(metadata)) {
    if (!all(sample_ids %in% rownames(metadata))) {
      stop("Sample names in expr_data must match rownames of metadata.")
    }

    metadata <- metadata[sample_ids, , drop = FALSE]

    anno_df <- data.frame(row.names = sample_ids)
    col_colors <- list()
    row_colors <- list()

    if (!is.null(ancestry_column)) {
      ancestry <- metadata[[ancestry_column]]
      ancestry_levels <- unique(ancestry)
      ancestry_colors <- setNames(c("navy", "goldenrod"), ancestry_levels)
      anno_df[[ancestry_column]] <- ancestry
      col_colors[[ancestry_column]] <- ancestry_colors
      row_colors[[ancestry_column]] <- ancestry_colors
    }

    if (!is.null(output_column)) {
      output <- metadata[[output_column]]
      output_levels <- unique(output)
      output_colors <- setNames(c("#027c58", "purple"), output_levels)
      anno_df[[output_column]] <- output
      col_colors[[output_column]] <- output_colors
      row_colors[[output_column]] <- output_colors
    }

    if (ncol(anno_df) > 0) {
      common_legend <- list(
        direction = "horizontal"
      )

      top_anno <- ComplexHeatmap::HeatmapAnnotation(
        df = anno_df,
        col = col_colors,
        show_annotation_name = FALSE,
        annotation_legend_param = common_legend
      )

      left_anno <- ComplexHeatmap::rowAnnotation(
        df = anno_df,
        col = row_colors,
        show_annotation_name = FALSE,
        annotation_legend_param = common_legend
      )
    }
  }

  heatmap_obj <- ComplexHeatmap::Heatmap(
    matrix            = cor_mat,
    name              = paste(method, "(r)"),
    col               = color_scale,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_cols,
    show_row_names    = show_names,
    show_column_names = show_names,
    row_names_gp      = grid::gpar(fontsize = point_size * 10),
    column_names_gp   = grid::gpar(fontsize = point_size * 10),
    top_annotation    = top_anno,
    left_annotation   = left_anno,
    heatmap_legend_param = list(
      direction = "horizontal"
    ),
    ...
  )

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

  ComplexHeatmap::draw(
    heatmap_obj,
    heatmap_legend_side     = "bottom",
    annotation_legend_side  = "bottom",
    show_heatmap_legend     = TRUE,
    show_annotation_legend  = TRUE,
    merge_legend            = TRUE
  )

  invisible(cor_mat)
}
