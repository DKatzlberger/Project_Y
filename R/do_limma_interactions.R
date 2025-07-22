do_limma_interactions <- function(
  Y,           # Expression matrix: samples × genes
  M,           # Sample metadata: samples × variables
  a_col,       # Ancestry variable
  g_col,       # Phenotype/condition variable
  a_train,     # Reference ancestry (e.g., "EA")
  a_infer,     # Comparison ancestry (e.g., "AA")
  g_levels,    # Phenotype levels: c(P₁, P₂)
  normalize = TRUE  # Whether to normalize expression via voom
) {
  clean <- function(x) gsub("-", "_", x)
  a_train <- clean(a_train)
  a_infer <- clean(a_infer)
  g_levels <- clean(g_levels)

  # Create group variable based on ancestry × phenotype
  M <- as.data.table(M)
  M[, group := factor(
    paste(clean(get(a_col)), clean(get(g_col)), sep = "."),
    levels = c(
      paste(a_train, g_levels[1], sep = "."),
      paste(a_train, g_levels[2], sep = "."),
      paste(a_infer, g_levels[1], sep = "."),
      paste(a_infer, g_levels[2], sep = ".")
    )
  )]

  # Design matrix
  X <- model.matrix(~ 0 + group, data = M)
  colnames(X) <- levels(M$group)

  # Normalize (if specified)
  V <- if (normalize) voom(t(Y), X) else list(E = t(Y), design = X)

  # Fit model
  fit_means <- lmFit(V$E, X)
  fit_means <- eBayes(fit_means)

  B <- colnames(X)

  # Algebraic contrasts
  C_math <- list(
    base_P1  = glue("{B[1]} - {B[3]}"),                         # β_Atrain_P1 - β_Ainfer_P1
    base_P2  = glue("{B[2]} - {B[4]}"),                         # β_Atrain_P2 - β_Ainfer_P2
    eff_Atr  = glue("{B[1]} - {B[2]}"),                         # β_Atrain_P1 - β_Atrain_P2
    eff_Ainf = glue("{B[3]} - {B[4]}"),                         # β_Ainfer_P1 - β_Ainfer_P2
    inter    = glue("({B[1]} - {B[2]}) - ({B[3]} - {B[4]})")    # interaction
  )

  # More math-like contrast labels
  C_labels <- c(
    glue("B_{a_train}_{g_levels[1]} - B_{a_infer}_{g_levels[1]}"),
    glue("B_{a_train}_{g_levels[2]} - B_{a_infer}_{g_levels[2]}"),
    glue("B_{a_train}_{g_levels[1]} - B_{a_train}_{g_levels[2]}"),
    glue("B_{a_infer}_{g_levels[1]} - B_{a_infer}_{g_levels[2]}"),
    glue("deltaB_{a_train} - deltaB_{a_infer}")  # interaction
  )

  # Build contrast matrix
  C <- makeContrasts(contrasts = C_math, levels = X)
  colnames(C) <- C_labels

  # Fit contrast model
  fit_contrasts <- contrasts.fit(fit_means, C)
  fit_contrasts <- eBayes(fit_contrasts)

  # Return
  list(
    fit_means     = fit_means,
    fit_contrasts = fit_contrasts,
    norm_matrix   = V$E
  )
}
