bootstrap_interaction <- function(
  X,          # Expression matrix 1
  Y,          # Expression matrix 2
  MX,         # Meta data 1
  MY,         # Meta data 2
  g_col,      # Condition column
  a_col,      # Ancestry column
  B = 1000,   # Number of iterations 
  seed        # Reproducibility 
){

  set.seed(seed)

  # Define groups that are compared
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]

  # Validate that the groups are factors with 2 identical levels
  if (!is.factor(g_X) || !is.factor(g_Y)) {
    stop("Group column must be a factor.")
  }

  if (!identical(levels(g_X), levels(g_Y))) {
    stop("Group levels must match between MX and MY.")
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

  if (length(a_X) != 1 || length(a_Y) != 1) {
    stop("Each metadata frame must contain only one ancestry level.")
  }

  # Difference in groups
  delta_X <- colMeans(X[g_X == g1, , drop = FALSE]) - colMeans(X[g_X == g2, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g1, , drop = FALSE]) - colMeans(Y[g_Y == g2, , drop = FALSE])

  # Interactions (difference of differences)
  T_obs <- delta_X - delta_Y

  # Combine the data
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)
  g_XY <- MXY[[g_col]]

  # Sample sizes
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Initilaize boot matrix
  T_boot <- matrix(0, nrow = B, ncol = ncol(XY))
  colnames(T_boot) <- colnames(XY)

  for (b in 1:B){
    # Randomly split into two pseudo ancestries
    indices <- sample(1:nrow(XY))
    p1 <- indices[1:n1]
    p2 <- indices[(n1 + 1):(n1 + n2)]

    # Group labels for each pseudo ancestry
    g1_p1 <- XY[p1[g_XY[p1] == g1], , drop = FALSE]
    g2_p1 <- XY[p1[g_XY[p1] == g2], , drop = FALSE]

    g1_p2 <- XY[p2[g_XY[p2] == g1], , drop = FALSE]
    g2_p2 <- XY[p2[g_XY[p2] == g2], , drop = FALSE]

    # Skip iteration if any group is empty
    if (any(sapply(list(g1_p1, g2_p1, g1_p2, g2_p2), nrow) == 0)) {
      T_boot[b, ] <- NA
      next
    }

    # Compute group effects and interaction
    d1 <- colMeans(g1_p1) - colMeans(g2_p1)
    d2 <- colMeans(g1_p2) - colMeans(g2_p2)
    T_boot[b, ] <- d1 - d2
  }

  # Remove any NAs caused by bad splits
  valid_rows <- complete.cases(T_boot)
  T_boot <- T_boot[valid_rows, , drop = FALSE]
  B_actual <- nrow(T_boot)

  # Empirical p-values
  p_vals <- sapply(1:ncol(XY), function(j) {
    mean(abs(T_boot[, j]) >= abs(T_obs[j]))
  })

  # Compute 95% CIs for each feature using percentiles
  ci_bounds <- apply(T_boot, 2, function(x) quantile(x, c(0.025, 0.975)))
  ci_lower <- ci_bounds[1, ]
  ci_upper <- ci_bounds[2, ]

  # FDR correction
  fdr <- p.adjust(p_vals, method = "BH")

  return(list(
    observed = data.table(
      feature = colnames(XY),
      T_obs = T_obs,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_vals,
      adj_p_value = fdr,
      coef = rep(glue("deltaB_{as.character(a_X)} - deltaB_{as.character(a_Y)}"), length(T_obs)),
      row.names = NULL
    ),
    T_boot = T_boot,
    B_used = B_actual
  ))
}
