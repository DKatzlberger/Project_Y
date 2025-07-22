bootstrap_correlation <- function(
  X,                  # Expression matrix 1 
  Y,                  # Expression matrix 2
  Z,                  # Expression matrix 3
  MX,                 # Meta data 1
  MY,                 # Meta data 2
  MZ,                 # Meta data 3
  g_col,              # Condition column
  a_col,              # Ancestry column
  method = "pearson", # Correaltion coefficient
  B = 1000,           # Number of iterations 
  seed                # Reproducibility 
){

  set.seed(seed)

  # Define groups that are compared
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]
  g_Z <- MZ[[g_col]]

  # Validate that the groups are factors with 2 identical levels
  if (!is.factor(g_X) || !is.factor(g_Y) || !is.factor(g_Z)) {
    stop("Group column must be a factor.")
  }

  if (!identical(levels(g_X), levels(g_Y)) || !identical(levels(g_X), levels(g_Z))) {
    stop("Group levels must match between MX, MY, and MZ.")
  }

  if (length(levels(g_X)) != 2) {
    stop("Group column must have exactly 2 levels.")
  }

  # Make sure the reference is correct
  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  # Validate ancestry is singular per dataset
  a_X <- unique(MX[[a_col]])
  a_Y <- unique(MY[[a_col]])
  a_Z <- unique(MZ[[a_col]])

  if (length(a_X) != 1 || length(a_Y) != 1 || length(a_Z) != 1) {
    stop("Each metadata frame must contain only one ancestry level.")
  }

  # Difference in groups
  delta_X <- colMeans(X[g_X == g1, , drop = FALSE]) - colMeans(X[g_X == g2, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g1, , drop = FALSE]) - colMeans(Y[g_Y == g2, , drop = FALSE])
  delta_Z <- colMeans(Z[g_Z == g1, , drop = FALSE]) - colMeans(Z[g_Z == g2, , drop = FALSE])

  # Correlation to Z
  cor_XZ <- cor(delta_X, delta_Z, use = "complete.obs", method = method)
  cor_YZ <- cor(delta_Y, delta_Z, use = "complete.obs", method = method)

  # Difference in correlation
  T_obs <- cor_XZ - cor_YZ

  # Bootstrap CI for individual correlations
  boot_cor_XZ <- replicate(B, {
    i <- sample(1:nrow(X), replace = TRUE)
    gX_i <- g_X[i]
    d <- colMeans(X[i[gX_i == g1], , drop = FALSE]) - colMeans(X[i[gX_i == g2], , drop = FALSE])
    cor(d, delta_Z, use = "complete.obs", method = method)
  })

  boot_cor_YZ <- replicate(B, {
    i <- sample(1:nrow(Y), replace = TRUE)
    gY_i <- g_Y[i]
    d <- colMeans(Y[i[gY_i == g1], , drop = FALSE]) - colMeans(Y[i[gY_i == g2], , drop = FALSE])
    cor(d, delta_Z, use = "complete.obs", method = method)
  })

  ci_XZ <- quantile(boot_cor_XZ, c(0.025, 0.975))
  ci_YZ <- quantile(boot_cor_YZ, c(0.025, 0.975))

  # Data frame for correlation results
  correlation_df <- data.table(
    method = method,
    target = c("X-Z", "Y-Z"),
    correlation = c(cor_XZ, cor_YZ),
    ci_lower = c(ci_XZ[1], ci_YZ[1]),
    ci_upper = c(ci_XZ[2], ci_YZ[2])
  )

  # Permutation test
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)
  g_XY <- MXY[[g_col]]

  n1 <- nrow(X)  # Size of ancestry X
  n2 <- nrow(Y)  # Size of ancestry Y

  # Initialize matrix for permutation results
  T_boot <- numeric(B)

  for (b in 1:B) {
    # Shuffle all sample indices
    indices <- sample(1:nrow(XY))
    p1 <- indices[1:n1]                # pseudo ancestry 1
    p2 <- indices[(n1 + 1):(n1 + n2)]  # pseudo ancestry 2

    # Get group labels for each pseudo ancestry
    g1_p1 <- XY[p1[g_XY[p1] == g1], , drop = FALSE]
    g2_p1 <- XY[p1[g_XY[p1] == g2], , drop = FALSE]

    g1_p2 <- XY[p2[g_XY[p2] == g1], , drop = FALSE]
    g2_p2 <- XY[p2[g_XY[p2] == g2], , drop = FALSE]

    # Skip if any group is empty
    if (any(sapply(list(g1_p1, g2_p1, g1_p2, g2_p2), nrow) == 0)) {
      T_boot[b] <- NA
      next
    }

    # Compute pseudo deltas
    d1 <- colMeans(g1_p1) - colMeans(g2_p1)
    d2 <- colMeans(g1_p2) - colMeans(g2_p2)

    # Compute correlation with reference delta
    cor1 <- cor(d1, delta_Z, use = "complete.obs", method = method)
    cor2 <- cor(d2, delta_Z, use = "complete.obs", method = method)

    # Store correlation difference
    T_boot[b] <- cor1 - cor2
  }

  # Clean up NAs
  T_boot <- T_boot[!is.na(T_boot)]
  B_actual <- length(T_boot)

  # Empirical pvalue
  p_val <- mean(abs(T_boot) >= abs(T_obs))
  fdr <- p.adjust(p_val, method = "BH")

  # CI for delta correlation difference (T_obs)
  ci_vals <- as.numeric(quantile(T_boot, c(0.025, 0.975)))
  ci_lower <- ci_vals[1]
  ci_upper <- ci_vals[2]

  observed <- data.table(
    feature = "global",  
    method = method,
    T_obs = T_obs,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_val,
    p_adj_value = fdr,
    contrast = paste0("Correlation delta with Z: ", a_X, " vs ", a_Y)
  )

  return(list(
    observed = observed,
    correlation = correlation_df,
    T_boot = T_boot,
    B_used = B_actual
  ))
}