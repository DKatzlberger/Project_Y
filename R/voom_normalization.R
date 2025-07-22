library(limma)
library(anndata)

voom_normalization <- function(adata) {
  # Apply voom normalization (voom expects features x samples)
  voom_result <- voom(t(adata$X))

  # Transpose result back to samples x features
  normalized_matrix <- t(voom_result$E)

  # Create new AnnData object with normalized expression
  new_adata <- AnnData(
    X = normalized_matrix,
    obs = adata$obs,
    var = adata$var
  )

  # Return a list with both old and new objects
  return(list(
    raw = adata,
    voom = new_adata
  ))
}
