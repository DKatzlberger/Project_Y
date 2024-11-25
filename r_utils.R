make.dir <- function(fp, overwrite = FALSE) {
  # Creates a directory 
  # can be specified if you want to overwrite.
  # fb: Path to directory.
  # overwrite: If the existing directory should be overwritten.
  if(!file.exists(fp)) {  
    make.dir(dirname(fp))
    dir.create(fp)
  } else {   if(isTRUE(overwrite)){
    unlink(fp, recursive = TRUE)
    dir.create(fp)
    print('Overwritting existing directory')
  } else {
    print('File exists!')
    } 
  }
} 

normalize_log <- function(X, e = 0.01){
  # Log normalization.
  # X: Vector or Matrix to be transformed.
  # e: threshold to correct for infinity values.
  X = log2(X + e)
}

create_design <- function(output_column, meta){
  # Creates design matrix 
  # Currently supports means model.
  # output_column: Condition to be calculated.
  
  # Design matrix
  design <- model.matrix(~0 + eval(as.name(output_column)), data = meta) 
  # Rename the columns of design matrix
  repl_str <- 'eval(as.name(output_column))'
  colnames(design) <- gsub(repl_str, '', colnames(design), fixed=TRUE)
  # Return
  design
}

create_contrast <- function(coef){
  # Creates contrast matrix.
  # Compares the average of one condition to the average of all other.
  # coef: Coefficients of the design matrix.
  
  # Contrast matrix
  coef1 <- coef
  contrast.matrix <- matrix(0, ncol = length(coef1), nrow=length(coef1))
  row.names(contrast.matrix) <- coef1
  colnames(contrast.matrix) <- coef1
  
  cx <- coef1[1]
  for(cx in coef1){
    colx <- cx
    # colx <- str_replace(cx, 'eval(as.name(setup$classification$output_column))', '')
    contrast.matrix[coef1, colx] <- 1
    contrast.matrix[row.names(contrast.matrix) != cx, colx] <- -1 * (1/(length(coef1)-1))
  }
  # Return
  contrast.matrix
}

extract_results <- function(limmaFit){
  # Extracts results from a topTable for all coefficients.
  # limmaFit: Results from limma fitted model.
  
  limmaRes <- list()
  for (x in colnames(coef(limmaFit))){
    # Extract for each coefficient all genes
    limmaRes[[x]] <- topTable(limmaFit, coef=x, number = Inf) |> 
      rownames_to_column('Feature')
  }
  limmaRes <- bind_rows(limmaRes, .id = "coef") 
  limmaRes <- filter(limmaRes, coef != "(Intercept)") 
  # Return
  limmaRes
}

visual_validation <- function(meta, 
  signal, 
  mean_stats, 
  contrast_stats, 
  goi,
  data_output_column){
  # To visually validate the output of stats. models.
  # meta: Meta data. 
  # signal: Signal of the molecular data.
  # mean_stats: Results from means model.
  # contrast_stats: Results from comparison.
  # goi: Genes to validate.

  mean_stats <- 
    mean_stats |> 
    as_tibble() |> 
    group_by(coef) |> 
    filter(Feature %in% goi) 

  contrast_stats <- 
    contrast_stats |> 
    as_tibble() |> 
    group_by(coef) |> 
    filter(Feature %in% goi) 
  
  signal <- 
    signal |> 
    as_tibble(rownames = 'Feature') |> 
    filter(Feature %in% goi) |> 
    as.data.frame() |> 
    column_to_rownames('Feature') |> 
    t() |> 
    as_tibble(rownames = 'sampleId') 
  
  mean_signal <- 
    meta |> 
    as_tibble(rownames = 'sampleId') |> 
    select(all_of(data_output_column), sampleId) |> 
    merge(signal, by = 'sampleId') |> 
    group_by_at(data_output_column) |> # interestng to use group_by_at
    summarise(across(where(is.numeric), mean)) |> 
    pivot_longer(cols = !all_of(data_output_column),
                 values_to = 'E',
                 names_to = 'Feature')
  
  # Common plot settings
  x_lab <- xlab('Comparison')
  y_lab <- ylab('Feature')

  # Plot stats. results
  contrast_stats_plot <- ggplot(
    contrast_stats,
    aes(
      x = coef,
      y = Feature,
      color = logFC,
      size = pmin(5, -log10(adj.P.Val))
    )
  ) +
       geom_point() +
       scale_color_gradient2(
        low="blue",
        high="red"
  ) +
  ggtitle('Results of comparison') +
  y_lab +
  x_lab +
  theme(axis.text.x = element_text(angle = 90))

    mean_stats_plot <- ggplot(
    mean_stats,
    aes(
      x = coef,
      y = Feature,
      color = logFC,
      size = pmin(5, -log10(adj.P.Val))
    )
  ) +
       geom_point() +
       scale_color_gradient2(
        low="blue",
        high="red"
  ) +
  ggtitle('Results of means model') +
  y_lab +
  x_lab +
  theme(axis.text.x = element_text(angle = 90))

  # Plot mean signal
  mean_signal_plot <- ggplot(
    mean_signal,
    aes(
      x = get(data_output_column),
      y = Feature,
      fill = E
    )
  ) +
  geom_tile() +
  scale_fill_gradient2(
        low="blue",
        high="red"
  ) +
  ggtitle('Mean signal in the data') +
  y_lab +
  x_lab +
  theme(axis.text.x = element_text(angle = 90))

  # Patchwork
  validation_plot <- contrast_stats_plot + mean_stats_plot + mean_signal_plot
  validation_plot
}

is_yml_file <- function(file_path) {
  # Check if the file exists and if the file has a .yml extension
  if (!file.exists(file_path)) {
    stop('Error: The file does not exist.')
  }
  if (tolower(tools::file_ext(file_path)) != 'yml') {
    stop('Error: The file is not a YML file.')
  }
}

compute_correlation <- function(data, method = "spearman") {
  # Compute correlation matrix
  correlation_matrix <- data |> 
    select(starts_with("logFC")) |> 
    cor(method = method)
  
  # Melt the upper triangle of the correlation matrix into a data.table
  melted_correlation <- as.data.table(as.table(correlation_matrix))[V1 < V2]
  
  # Rename the correlation value column
  setnames(melted_correlation, old = "N", new = str_to_title(method))
  
  return(head(melted_correlation, 2))
}
