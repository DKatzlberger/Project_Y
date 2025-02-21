variance_explained <- data.table()
for (condition in unique(adata$obs[[output_column]])){
  # Meta data
  filtered_meta <- adata$obs |>
    as_tibble(rownames = "Idx") |>
    filter(!!sym(output_column == condition)) |>
    dplyr::select(Idx, all_of(ancestry_column))

  # Expression
  filtered_expression <- adata$X |>
    as_tibble(rownames = "Idx") |>
    inner_join(filtered_meta, by = "Idx")

  # Reshape
  filtered_expression |>
    pivot_longer(
      -c(Idx, all_of(ancestry_column)),
      names_to = "Feature",
      values_to = "Expression"
    )
  
  # R2 per condition
  con_variance_explained <- calculate_variance_explained(
    expression_df = filtered_expression,
    check_varaince = c(ancestry_column),
    batch_size = 1000
  ) |>
  mutate(Condition = condition)

  # Combine
  variance_explained <- bind_rows(variance_explained, con_variance_explained)
}



variance_across_ancestry <- data.table()
for (condition in unique(adata$obs[[output_column]])){
  # Meta data
  filtered_meta <- adata$obs |>
    as_tibble(rownames = "Idx") |>
    filter(!!sym(output_column == condition)) |>
    dplyr::select(Idx, all_of(ancestry_column))
  
  # Expression
  filtered_expression <- adata$X |>
    as_tibble(rownames = "Idx") |>
    inner_join(filtered_meta, by = "Idx")

  # Calculate variance
  con_variance_across_ancestry <- filtered_expression |>
    group_by(!!sym(ancestry_column), Feature) |>
    summarize(
      Variance = var(Expression, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(Condition == condition)
  
  # Combine
  variance_across_ancestry <- bind_rows(variance_across_ancestry, con_variance_across_ancestry)
}



residual_threshold <- 1
upper_threshold <- (intercept + slope * observed_values + residual_threshold)[, 1]
lower_threshold <- (intercept + slope * observed_values - residual_threshold)[, 1]

scatter_plot_data <- baseline |>
  mutate(Comparison = paste(Comparison, Condition)) |>
  dplyr::select(Feature, Comparison, logFC) |>
  pivot_wider(names_from = Comparison, values_from = logFC) |>
  mutate(
    abs_residual = abs_residuals,
    above_threshold = abs_residual > residual_threshold,
    upper_threshold = intercept + slope * observed_values + residual_threshold,
    lower_threshold = intercept + slope * observed_values - residual_threshold
  )