do_interactions <- function(
  expr_data,
  sample_metadata,
  ancestry_column,
  output_column,
  train_ancestry,
  infer_ancestry,
  comparison
) {
  library(limma)
  library(data.table)
  library(glue)

  # Replace "-" with "_" for valid column names
  clean <- function(x) gsub("-", "_", x)
  train_ancestry <- clean(train_ancestry)
  infer_ancestry <- clean(infer_ancestry)
  comparison <- clean(comparison)

  # Format metadata and build design matrix
  sample_metadata <- as.data.table(sample_metadata)
  sample_metadata[
    ,
    group := factor(
      paste(clean(get(ancestry_column)), clean(get(output_column)), sep = "."),
      levels = c(
        paste(train_ancestry, comparison[1], sep = "."),
        paste(train_ancestry, comparison[2], sep = "."),
        paste(infer_ancestry, comparison[1], sep = "."),
        paste(infer_ancestry, comparison[2], sep = ".")
      )
    )
  ]

  design <- model.matrix(~ 0 + group, data = sample_metadata)
  colnames(design) <- levels(sample_metadata$group)

  # Normalize expression data
  v <- voom(expr_data, design)

  # Fit means model
  fit_means <- lmFit(v, design)
  fit_means <- eBayes(fit_means)

  # Define human-readable contrast names
  contrast_terms <- list(
    baseline_1      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}.Baseline"),
    baseline_2      = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[2]}.Baseline"),
    relationship_1  = glue("{train_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relationship"),
    relationship_2  = glue("{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Relationship"),
    interaction     = glue("{train_ancestry}_vs_{infer_ancestry}.{comparison[1]}_vs_{comparison[2]}.Interaction")
  )

  # Use design colnames for contrast math
  cols <- colnames(design)
  contrast_calculations <- list(
    baseline_1      = glue("{cols[1]} - {cols[3]}"),
    baseline_2      = glue("{cols[2]} - {cols[4]}"),
    relationship_1  = glue("{cols[2]} - {cols[1]}"),
    relationship_2  = glue("{cols[4]} - {cols[3]}"),
    interaction     = glue("({cols[2]} - {cols[1]}) - ({cols[4]} - {cols[3]})")
  )

  # Build contrast matrix
  contrast_matrix <- makeContrasts(
    contrasts = contrast_calculations,
    levels    = design
  )

  colnames(contrast_matrix) <- contrast_terms

  # Fit contrasts
  fit_contrasts <- contrasts.fit(fit_means, contrast_matrix)
  fit_contrasts <- eBayes(fit_contrasts)

  # Return both fits
  list(
    fit_means      = fit_means,
    fit_contrasts  = fit_contrasts
  )
}
