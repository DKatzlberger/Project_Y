is_yaml_file <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("Error: The file does not exist.")
  }

  # Check if the file has a valid YAML extension (.yml or .yaml)
  ext <- tolower(tools::file_ext(file_path))
  if (!(ext %in% c("yml", "yaml"))) {
    stop("Error: The file is not a YAML file. Please ensure the file has a .yml or .yaml extension.")
  }
}

# Normalization functions:
normalize_log <- function(X, e = 1, ...){
  # Log normalization.
  # X: Vector or Matrix to be transformed.
  # e: threshold to correct for infinity values.
  return(log2(X + e))
}

beta_to_mvalue <- function(betas, epsilon = 0.00001, ...) {
  # This function transforms beta values (ranging from 0 to 1) into M-values 
  # using a logit transformation. M-values provide a more statistically normal 
  # distribution compared to beta values and are commonly used in DNA methylation 
  # analysis.

  # betas: A numeric vector of beta values (ranging between 0 and 1).
  # epsilon: A small numeric value (default: 0.00001) to prevent extreme 
  # log transformations by clamping values close to 0 or 1.
	if (!is.numeric(betas)) {
		stop("invalid value for betas")
	}
	if (!(is.numeric(epsilon) && length(epsilon) == 1 && (!is.na(epsilon)))) {
		stop("invalid value for epsilon")
	}
	if (epsilon < 0 || epsilon > 0.5) {
		stop("invalid value for epsilon; expected 0 <= epsilon <= 0.5")
	}
	betas[betas < epsilon] <- epsilon
	betas[betas > (1 - epsilon)] <- 1 - epsilon
	return(log2(betas / (1 - betas)))
}

voom_normalization <- function(count_matrix, design_matrix) {
  # Ensure input is a matrix
  if (!is.matrix(count_matrix)) {
    stop("Error: Input count_matrix must be a matrix.")
  }
  
  # Create a DGEList object
  dge <- DGEList(counts = count_matrix)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge)
  
  # Apply voom transformation
  voom_res <- voom(dge, design_matrix, plot = FALSE)
  
  return(voom_res)
}

# Normalization lookup list
normalization_methods <- list(
  "transcriptomics" = list(
    "limma_voom" = list(
      "function" = function(X, design, ...) voom_normalization(X, design, ...),  
      "output_name" = "logCPM (voom)"
    ),
    "normalize_log" = list(
      "function" = function(X, ...) normalize_log(X),
      "output_name" = "log2-transformed values"
    ),
    "normalize_zscore" = list(
      "function" = function(X, ...) scale(X),
      "output_name" = "Z-score normalized values"
    ),
    "raw" = list(
      "function" = function(X, ...) X,
      "output_name" = "Input values"
    )
  ),
  "methylation" = list(
    "beta_to_mvalue" = list(
      "function" = function(X, ...) beta_to_mvalue(X),
      "output_name" = "M-values"
    ),
    "normalize_log" = list(
      "function" = function(X, ...) normalize_log(X),
      "output_name" = "log2-transformed values"
    ),
    "normalize_zscore" = list(
      "function" = function(X, ...) scale(X),
      "output_name" = "Z-score normalized values"
    ),
    "raw" = list(
      "function" = function(X, ...) X,
      "output_name" = "Input values"
    )
  ),
  "proteomics" = list(
    "normalize_zscore" = list(
      "function" = function(X, ...) scale(X),
      "output_name" = "Z-score normalized values"
    ),
    "raw" = list(
      "function" = function(X, ...) X,
      "output_name" = "Input values"
    ),
    "normalize_log" = list(
      "function" = function(X, ...) normalize_log(X),
      "output_name" = "log2-transformed values"
    )
  )
)




# Filter features
# TMM Norm factors
calculate_tmm_norm_factors <- function(data) {
  # 'data' is a matrix with samples as rows and genes as columns (samples x genes)
  
  # Transpose the data to get genes as rows and samples as columns (for easier calculation)
  data_transposed <- t(data)
  
  # Step 1: Compute the geometric mean of the counts for each gene (row)
  geometric_mean <- apply(data_transposed, 1, function(x) exp(mean(log(x + 1))))  # log-transformed for stability
  
  # Step 2: Compute M-values (log-ratio) between each sample and the geometric mean
  # M-value = log2(count_in_sample / geometric_mean)
  M_values <- sweep(data_transposed, 1, geometric_mean, FUN = "/")  # Element-wise division
  M_values <- log2(M_values + 1)  # Add pseudocount for stability
  
  # Step 3: Trim the extreme M-values to reduce the influence of highly differentially expressed genes
  # Here we trim the top and bottom 5% of M-values
  trim_percent <- 0.05
  M_values_trimmed <- apply(M_values, 2, function(x) {
    lower <- quantile(x, trim_percent)
    upper <- quantile(x, 1 - trim_percent)
    x[x < lower] <- lower
    x[x > upper] <- upper
    return(x)
  })
  
  # Step 4: Calculate the normalization factors
  # For each sample, calculate the median of the trimmed M-values for all genes
  norm_factors <- apply(M_values_trimmed, 2, median)
  
  # Normalize by the median normalization factor
  norm_factors <- norm_factors / median(norm_factors)
  
  return(norm_factors)
}

# CPM
cpm <- function(data, norm_factors = NULL, log = FALSE) {
  # Converts raw counts to CPM (Counts Per Million), optionally using normalization factors.
  # Optionally, log-transform the CPM values.
  # Args:
  #   data (matrix): Raw counts (samples x genes).
  #   norm_factors (vector, optional): Normalization factors for each sample. Default is NULL.
  #   log (bool): Whether to log-transform the CPM values (default is FALSE).
  # Returns:
  #   matrix: CPM or logCPM values (samples x genes).
  
  # Check if the number of samples in the data matches the length of norm_factors
  if (!is.null(norm_factors) && nrow(data) != length(norm_factors)) {
    stop("The number of samples must match the length of the normalization factors.")
  }
  
  # If norm_factors are provided, apply them to the raw counts
  if (!is.null(norm_factors)) {
    # Apply normalization factors to adjust raw counts (divide by norm factors)
    data <- sweep(data, 1, norm_factors, FUN = "/")
  }
  
  # Calculate total reads for each sample (sum of all gene counts in each sample)
  total_reads_per_sample <- rowSums(data)
  
  # Convert raw counts to CPM
  cpm_data <- sweep(data, 1, total_reads_per_sample, FUN = "/") * 1e6
  
  # Optionally, log-transform the CPM data
  if (log) {
    cpm_data <- log2(cpm_data + 1)  # Add pseudocount to avoid log(0)
  }
  
  return(cpm_data)
}

# By signal
signal_by_percentile <- function(data, percentile) {
  # Calculates the signal threshold for genes based on the specified percentile of their total counts.
  # Args:
  #   data (matrix): Expression data (samples x genes).
  #   percentile (numeric): Percentile to calculate the signal threshold (e.g., 25 for 25th percentile).
  # Returns:
  #   numeric: The calculated signal threshold for genes based on the specified percentile.
  
  # Calculate total signal for each gene across all samples (sum along rows)
  gene_signal <- colSums(data)
  
  # Calculate the threshold using the specified percentile
  signal <- quantile(gene_signal, probs = percentile / 100)
  
  return(signal)
}

filter_by_signal <- function(data, min_signal = 10, max_signal = NULL, min_samples_ratio = 0.5) {
  # Filters genes based on minimum and maximum signal and presence in at least a given fraction of samples.
  # Args:
  #   data (matrix): Expression data (samples x genes).
  #   min_signal (numeric): Minimum total signal required for a gene (default: 10).
  #   max_signal (numeric): Maximum total signal for a gene (optional).
  #   min_samples_ratio (numeric): Minimum fraction of samples a gene must be expressed in (default: 50%).
  # Returns:
  #   logical: A logical vector indicating which genes were retained.
  
  # Calculate the number of samples
  num_samples <- nrow(data)
  
  # Calculate the minimum number of samples for the given ratio
  min_samples <- floor(min_samples_ratio * num_samples)
  
  # Calculate total signal for each gene across all samples (sum along rows)
  gene_signal <- colSums(data)
  
  # Calculate how many samples express each gene (non-zero signal counts)
  gene_samples <- colSums(data > 0)
  
  # Filter genes based on min_signal, max_signal, and min_samples_ratio
  if (is.null(max_signal)) {
    # If max_signal is not provided, only filter by min_signal and sample count
    filtered_genes <- (gene_signal >= min_signal) & (gene_samples >= min_samples)
  } else {
    # If max_signal is provided, filter by both min_signal and max_signal
    filtered_genes <- (gene_signal >= min_signal) & (gene_signal <= max_signal) & (gene_samples >= min_samples)
  }
  
  return(filtered_genes)
}

# By variance
variance_by_percentile <- function(data, percentile = 25) {
  # Calculates the variance threshold for methylation genes based on the specified percentile.
  # Args:
  #   data (matrix): Methylation beta values (samples x genes).
  #   percentile (numeric): Percentile to calculate the variance threshold (default: 25).
  # Returns:
  #   numeric: The calculated variance threshold based on the specified percentile.
  
  # Compute variance for each gene across samples (variance along columns)
  gene_variances <- apply(data, 2, var)
  
  # Determine the variance threshold at the given percentile
  variance <- quantile(gene_variances, probs = percentile / 100)
  
  return(variance)
}

filter_by_variance <- function(data, var_threshold = 0.01, min_samples_ratio = 0.5) {
  # Filters methylation genes based on variance and presence in at least a given fraction of samples.
  # Args:
  #   data (matrix): Methylation beta values (samples x genes).
  #   var_threshold (numeric): Minimum variance threshold to retain a gene (default: 0.01).
  #   min_samples_ratio (numeric): Minimum fraction of samples where the gene must be methylated above 0 (default: 50%).
  # Returns:
  #   logical: A logical vector indicating which genes were retained.
  
  # Calculate the number of samples
  num_samples <- nrow(data)
  
  # Calculate the minimum number of samples for the given ratio
  min_samples <- floor(min_samples_ratio * num_samples)
  
  # Calculate variance for each gene across samples (variance along columns)
  gene_variances <- apply(data, 2, var)
  
  # Count how many samples have a methylation value > 0 for each gene (non-zero beta value)
  methylated_samples <- colSums(data > 0)
  
  # Filter genes based on variance threshold and the number of methylated samples
  filtered_genes <- (gene_variances > var_threshold) & (methylated_samples >= min_samples)
  
  return(filtered_genes)
}





# Function to check if covariate is present in settings
get_covariate <- function(setup) {
  # Function to extract the covariate if present and valid
  if ("covariate" %in% names(setup) && 
      !is.null(setup$covariate) &&
      any(nzchar(setup$covariate))) {
    return(setup$covariate)
  } else {
    return(NULL)
  }
}

classify_covariates <- function(df, covariates) {
  # Initialize lists to store continuous and discrete covariates
  classified_covariates <- list(
    continuous = c(),
    discrete = c()
  )
  
  # Iterate through each covariate
  for (covariate in covariates) {
    # Skip covariates not present in the DataFrame
    if (!covariate %in% colnames(df)) {
      next
    }
    
    # Extract the column data
    col_data <- df[[covariate]]
    unique_values <- length(unique(col_data))
    
    # Check if column is a factor or character (discrete)
    if (is.factor(col_data) || is.character(col_data)) {
      classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
    } else if (is.numeric(col_data)) {
      # Heuristic: Numerical columns with <= 10 unique values are discrete
      if (unique_values <= 10) {
        classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
      } else {
        classified_covariates$continuous <- c(classified_covariates$continuous, covariate)
      }
    } else {
      # Treat other data types as discrete by default
      classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
    }
  }
  
  return(classified_covariates)
}


check_covariate_conditions <- function(data, covariates) {
  
  # Iterate through the list of covariates
  for (covariate in covariates) {
    
    # Count occurrences of each unique value in the covariate
    counts <- data |> group_by_at(covariate) |> count()
    
    # Get the unique values directly from the dataset
    possible_values <- unique(data[[covariate]])
    
    # Check if there are at least two unique levels in the covariate
    if (length(possible_values) < 2) {
      stop(paste("Covariate", covariate, "has fewer than two unique values"))
    }
    
    # Check if all 'possible_values' are represented in the counts
    missing_values <- setdiff(possible_values, counts[[covariate]])
    if (length(missing_values) > 0) {
      stop(paste(
        "Not enough samples for:",
        paste(missing_values, collapse = ", "),
        "to calculate covariate for covariate",
        covariate
      ))
    }

    # Check if there are at least 2 samples per value
    # if (any(counts$n <= 2)) {
    #   stop(paste("Not enough samples per covariate value for:", covariate))
    # }
  }
}


# Define the function
create_stratification_plot <- function(datasets, to_visualize_columns) {
  # Combine datasets with their labels
  combined_data <- bind_rows(
    lapply(names(datasets), function(set_name) {
      datasets[[set_name]] |> mutate(Set = set_name)
    })
  )
  
  # Base plot for the first column in to_visualize_columns
  base_plot <- combined_data |>
    ggplot(aes(x = Set, fill = !!sym(to_visualize_columns[1]))) +
    geom_bar(position = "fill") +
    ylab("Proportion") +
    ggtitle(paste("Distribution of", to_visualize_columns[1]))
  
  # If there are additional columns, create plots for them
  if (length(to_visualize_columns) > 1) {
    covariate_plots <- lapply(to_visualize_columns[-1], function(column) {
      combined_data |>
        ggplot(aes(x = Set, fill = !!sym(column))) +
        geom_bar(position = "fill") +
        ylab("Proportion") +
        ggtitle(paste("Distribution of", column))
    })
    
    # Combine all plots using patchwork
    final_plot <- Reduce(`/`, c(list(base_plot), covariate_plots))
  } else {
    # Only the base plot if no additional columns are specified
    final_plot <- base_plot
  }
  
  # Return the final patchwork plot
  return(final_plot)
}


# Function to create design matrix
create_design <- function(output_column, meta, covariates = NULL) {
  # Creates design matrix
  # Supports means model with optional covariates.
  # output_column: Condition to be calculated.
  # covariates: Optional list of covariates to include.
  
  # Ensure all factor levels in `meta` are relevant (drop unused levels)
  meta <- meta |> mutate(across(where(is.factor), droplevels))
  
  if (!is.null(covariates)) {
    # Check that all covariates exist in the metadata
    missing_covariates <- setdiff(covariates, colnames(meta))
    if (length(missing_covariates) > 0) {
      stop(paste("Covariates", paste(missing_covariates, collapse = ", "), "are not present in the metadata."))
    }
    
    # Create the formula for the model matrix with multiple covariates
    formula <- as.formula(paste("~0 +", output_column, "+", paste(covariates, collapse = " + ")))
    design <- model.matrix(formula, data = meta)
  } else {
    # Design matrix without covariates
    formula <- as.formula(paste("~0 +", output_column))
    design <- model.matrix(formula, data = meta)
  }
  
  # Safely rename the columns of the design matrix
  # Remove only prefixes for `output_column` and `covariates`, while leaving the rest intact
  colnames(design) <- gsub(paste0("^", output_column, ""), "", colnames(design)) # Remove `output_column` prefix
  
  if (!is.null(covariates)) {
    for (cov in covariates) {
      if (is.factor(meta[[cov]])) {
        # Remove factor variable name from column names (e.g., "covLevel1" -> "Level1")
        colnames(design) <- gsub(paste0("^", cov, ""), "", colnames(design))
      }
    }
  }
  
  # Return the design matrix
  return(design)
}

# Create desing matrix with covariate if available
# create_design <- function(output_column, meta, covariate = NULL) {
#   # Creates design matrix
#   # Supports means model with optional covariate.
#   # output_column: Condition to be calculated.
#   # covariate: Optional covariate to include.
  
#   if (!is.null(covariate)) {
#     if (!(covariate %in% colnames(meta))) {
#       stop(paste("Covariate", covariate, "is not present in the metadata."))
#     }
#     # Design matrix including covariate
#     design <- model.matrix(~0 + eval(as.name(output_column)) + eval(as.name(covariate)), data = meta)
#   } else {
#     # Design matrix without covariate
#     design <- model.matrix(~0 + eval(as.name(output_column)), data = meta)
#   }
  
#   # Rename the columns of the design matrix
#   if (!is.null(covariate)) {
#     # Replace both output_column and covariate from column names
#     output_repl_str <- 'eval(as.name(output_column))'
#     covariate_repl_str <- 'eval(as.name(covariate))'
#     colnames(design) <- gsub(output_repl_str, '', colnames(design), fixed = TRUE)
#     colnames(design) <- gsub(covariate_repl_str, '', colnames(design), fixed = TRUE)
#   } else {
#     # Replace only output_column from column names
#     output_repl_str <- 'eval(as.name(output_column))'
#     colnames(design) <- gsub(output_repl_str, '', colnames(design), fixed = TRUE)
#   }
  
#   # Return the design matrix
#   return(design)
# }

create_contrast <- function(coef, conditions){
  # Creates contrast matrix.
  # coef: Coefficients of the design matrix.
  # conditions: A vector of coefficient names that represent conditions of interest.
  
  # Validate that conditions are a subset of coef
  if (!all(conditions %in% coef)) {
    stop("All conditions must be present in the coef vector.")
  }
  
  # Initialize the contrast matrix with zeros, with rows for all coefficients
  contrast.matrix <- matrix(0, ncol = length(coef), nrow = length(coef))
  
  # Assign proper row and column names to the contrast matrix
  row.names(contrast.matrix) <- coef
  colnames(contrast.matrix) <- coef
  
  # Loop through each condition of interest and create the contrasts
  for (cx in conditions) {
    contrast.matrix[cx, cx] <- 1  # Set the contrast for the condition to 1
    other_conditions <- setdiff(conditions, cx)  # Other conditions get -1
    contrast.matrix[other_conditions, cx] <- -1 * (1 / length(other_conditions))
  }
  
  # Ensure that all other coefficients not in conditions have their rows set to 0
  remaining_coeffs <- setdiff(coef, conditions)
  contrast.matrix[remaining_coeffs, ] <- 0

  # Remove all columns with only zero values
  contrast.matrix <- contrast.matrix[, colSums(contrast.matrix != 0) > 0]
  
  # Return the contrast matrix
  return(contrast.matrix)
}

# create_contrast <- function(coef){
#   # Creates contrast matrix.
#   # Compares the average of one condition to the average of all other.
#   # coef: Coefficients of the design matrix.
  
#   # Contrast matrix
#   coef1 <- coef
#   contrast.matrix <- matrix(0, ncol = length(coef1), nrow=length(coef1))
#   row.names(contrast.matrix) <- coef1
#   colnames(contrast.matrix) <- coef1
  
#   cx <- coef1[1]
#   for(cx in coef1){
#     colx <- cx
#     # colx <- str_replace(cx, 'eval(as.name(setup$classification$output_column))', '')
#     contrast.matrix[coef1, colx] <- 1
#     contrast.matrix[row.names(contrast.matrix) != cx, colx] <- -1 * (1/(length(coef1)-1))
#   }
#   # Return
#   contrast.matrix
# }


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


compute_correlation <- function(data, method = "spearman") {
  # Compute correlation matrix from DGE results
  correlation_matrix <- data |> 
    select(starts_with("logFC")) |> 
    cor(method = method)
  
  # Melt the upper triangle of the correlation matrix into a data.table
  melted_correlation <- as.data.table(as.table(correlation_matrix))[V1 < V2]
  
  # Rename the correlation value column
  setnames(melted_correlation, old = "N", new = str_to_title(method))
  
  return(melted_correlation)
}


# Function to calculate the mean response per output_column and Feature
calculate_observed_mean_response <- function(data, meta, output_column) {
  
  # Combine the data and metadata
  combined_data <- data |>
    left_join(meta, by = "patient_id") |>
    pivot_longer(cols = -c(all_of(output_column), patient_id), 
                 names_to = "Feature", 
                 values_to = "Value")
  
  # Calculate mean responses per feature and cancer type
  mean_responses <- combined_data |>
    group_by(!!sym(output_column), Feature) |>
    summarise(mean_response = mean(Value), .groups = "drop")
  
  return(mean_responses)
}


# Function to calculate contrast for each feature and format output
calculate_observed_contrast <- function(mean_responses, output_column, cancer_types = NULL) {
  
  # Calculate the mean response for each cancer_type and Feature
  mean_responses_grouped <- mean_responses |>
    group_by(Feature, !!sym(output_column)) |>
    summarise(mean_response = mean(mean_response), .groups = "drop")
  
  # If no specific cancer types are provided, use the unique cancer types from the data
  if (is.null(cancer_types)) {
    cancer_types <- unique(mean_responses[[output_column]])
  }
  
  # Initialize a list to store the contrast results
  contrast_results <- list()
  
  # Loop over all features to calculate contrasts
  for (feature in unique(mean_responses$Feature)) {
    
    # Subset the data for the specific feature
    feature_data <- mean_responses_grouped |> filter(Feature == feature)
    
    # Calculate the contrast for each pair of cancer types
    for (i in 1:(length(cancer_types) - 1)) {
      for (j in (i + 1):length(cancer_types)) {
        type1 <- cancer_types[i]
        type2 <- cancer_types[j]
        
        # Get the mean response for each cancer type
        mean1 <- feature_data$mean_response[feature_data[[output_column]] == type1]
        mean2 <- feature_data$mean_response[feature_data[[output_column]] == type2]
        
        # Compute the contrast (e.g., Type1 - Type2)
        contrast_1 <- mean1 - mean2
        contrast_2 <- mean2 - mean1  # Inverse comparison
        
        # Store the result in a tibble with comparison, feature, and contrast
        contrast_results[[paste(feature, type1, "vs", type2)]] <- tibble(
          Comparison = paste(type1),
          Feature = feature,
          Observed_contrast = contrast_1
        )
        
        contrast_results[[paste(feature, type2, "vs", type1)]] <- tibble(
          Comparison = paste(type2),
          Feature = feature,
          Observed_contrast = contrast_2
        )
      }
    }
  }
  
  # Combine all results into one data frame
  contrast_df <- bind_rows(contrast_results) |>
    arrange(Comparison)
  
  # Return the contrast data frame
  return(contrast_df)
}


# Interactions check variance
# Optimized function
calculate_variance_explained <- function(expression_df, check_variance, batch_size = 100, num_cores = detectCores() - 1) {
  setDT(expression_df)
  
  # Split the data by feature (gene) into batches
  feature_list <- unique(expression_df$Feature)
  batches <- split(feature_list, ceiling(seq_along(feature_list) / batch_size))
  
  # Set up parallel cluster
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("expression_df", "check_variance"), envir = environment())
  clusterEvalQ(cl, library(data.table))
  
  # Parallelize the processing of each batch
  batch_results <- parLapply(cl, batches, function(batch) {
    # Subset data for the current batch
    batch_data <- expression_df[Feature %in% batch, ]
    
    # Function to calculate R2 for each Feature and check_variance variable
    calculate_r2_for_feature <- function(.x, var) {
      var_data <- .x[[var]]
      
      # Determine if the variable is numeric or categorical
      if (is.numeric(var_data)) {
        # Numeric predictor (e.g., age)
        formula <- reformulate(var, response = "Expression")
        model <- lm(formula, data = .x)
        r2 <- summary(model)$r.squared
      } else {
        # Categorical predictor: Convert character to factor
        var_data <- as.factor(var_data)
        
        # Ensure at least two unique levels to avoid contrast error
        if (nlevels(var_data) > 1) {
          .x[[var]] <- var_data  # Assign back converted factor
          formula <- reformulate(var, response = "Expression")
          model <- lm(formula, data = .x)
          r2 <- summary(model)$r.squared
        } else {
          r2 <- NA  # Skip models with only one category
        }
      }
      return(data.table(Feature = unique(.x$Feature), Variable = var, R2 = r2))
    }
    
    # Efficiently calculate R2 for each variable in check_variance
    feature_results <- lapply(batch, function(feature) {
      feature_data <- batch_data[Feature == feature, ]
      
      # Calculate R2 for each variable in check_variance
      result <- lapply(check_variance, function(var) {
        calculate_r2_for_feature(feature_data, var)
      })
      
      # Combine the results for this feature
      rbindlist(result)
    })
    
    # Combine all features' results for this batch
    rbindlist(feature_results)
  })
  
  # Combine all batches' results efficiently
  variance_df <- rbindlist(batch_results)
  
  # Stop the cluster
  stopCluster(cl)
  
  return(variance_df)
}


# Permutation test
permutation_test <- function(data, value, group, paired = NULL, tail = "both") {
  
  # # Check if any group has identical values (zero variance)
  # var_check <- data |>
  #   group_by(!!sym(group)) |>
  #   summarise(var_value = var(!!sym(value)), .groups = 'drop')
  
  # # If any group has zero variance, return p-value of 1
  # if (any(var_check$var_value == 0)) {
  #   message("Perm test: One or more groups have identical values. Returning `NA`.")
  #   return(NA)  # Return the p-value of NA directly
  # }

  # Quote the column names for later use in evaluation
  value <- sym(value)  
  group <- sym(group)
  
  # Prepare the data by selecting relevant columns
  # If paired is provided, also select the paired column
  if (!is.null(paired)) {
    perm_data <- data |>
      select(!!value, !!group, !!sym(paired)) |>
      mutate(
        !!group := as.factor(!!group),
        !!sym(paired) := as.factor(!!sym(paired))  # Ensure the paired column is a factor
      )
  } else {
    perm_data <- data |>
      select(!!value, !!group) |>
      mutate(
        !!group := as.factor(!!group)
      )
  }
  
  # If paired is not NULL, use it as the paired column (blocking variable)
  if (!is.null(paired)) {
    # Construct the formula for paired test using the '|' operator for blocking
    formula <- as.formula(paste(deparse(value), "~", deparse(group), "|", paired))
  } else {
    # Construct the formula for the independence test (no pairing)
    formula <- as.formula(paste(deparse(value), "~", deparse(group)))
  }
  
  # Perform the permutation test (independence test)
  perm_test_result <- independence_test(formula = formula, data = perm_data)
  
  # Extract test statistic (Z-score)
  test_stat <- statistic(perm_test_result, type = "standardized")
  
  # Compute p-values
  if (tail == "right") {
    p_value <- as.numeric(1 - pnorm(test_stat))  # Right-tailed test
  } else if (tail == "left") {
    p_value <- as.numeric(pnorm(test_stat))  # Left-tailed test
  } else {
    p_value <- as.numeric(pvalue(perm_test_result))  # Default two-tailed test
  }
  
  return(p_value)
}


# permutation_test <- function(data, value_col, group_col) {
  
#   # Check if any group has identical values (zero variance)
#   var_check <- data |>
#   group_by(!!sym(group_col)) |>
#   summarise(var_value = var(!!sym(value_col)), .groups = 'drop')
  
#   # If any group has zero variance, return p-value of 1
#   if (any(var_check$var_value == 0)) {
#     message("Perm test: One or more groups have identical values. Returning `NA`.")
#     return(NA)  # Return the p-value of NA directly
#   }

#   # Quote the column names for later use in evaluation
#   value_col <- sym(value_col)  # Quoting the value column
#   group_col <- sym(group_col)  # Quoting the group column
  
#   # Prepare the data by selecting relevant columns and ensuring the group variable is a factor
#   perm_data <- data |>
#     select(!!value_col, !!group_col) |>
#     mutate(!!group_col := as.factor(!!group_col))  # Ensure group_col is a factor
  
#   # Construct the formula for the independence test
#   formula <- as.formula(paste(deparse(value_col), "~", deparse(group_col)))
  
#   # Perform the permutation test (independence test)
#   perm_test_result <- independence_test(formula = formula, data = perm_data)
  
#   # Return the result of the permutation test
#   p_value = as.numeric(pvalue(perm_test_result))
#   return(p_value)
# }

# Take a non stratified sample from a population by proportortions
# sample_with_proportions <- function(data, proportions, sizes, seed) {
#   # Set seed for reproducibility
#   set.seed(seed)
  
#   # List to store the sampled subsets
#   samples <- list()
  
#   # Sample data with replacement for each proportion
#   for (i in seq_along(sizes)) {
#     sample_name <- paste0("proportion_", gsub("\\.", "_", as.character(proportions[i])))
#     indexes <- sample(nrow(data), size = sizes[i], replace = FALSE)
#     samples[[sample_name]] <- data[indexes, ]
#   }
  
#   return(samples)
# }


# Define a generalized function to load and bind data from CSV files
fload_data <- function(folders, file_name) {
  combined_data <- data.table()  # Initialize an empty data.table

  # Loop through each folder to read and bind the specified CSV file
  for (folder in folders) {
    file <- file.path(folder, file_name)

    if (file.exists(file)) {
      # Read the data from the CSV file and bind it to the result
      data <- fread(file)
      combined_data <- rbindlist(list(combined_data, data), use.names = TRUE, fill = TRUE)
    } else {
      warning("File does not exist: ", file)
    }
  }

  # Return the combined data.table
  return(combined_data)
}




sample_by_size <- function(data, sizes, seed, class_column) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # List to store the sampled subsets
  samples <- list()
  
  # Ensure each class has at least two samples
  for (i in seq_along(sizes)) {
    sample_name <- paste0("size_", gsub("\\.", "_", as.character(sizes[i])))
    size <- sizes[i]
    
    # Split data by class
    class_groups <- split(data$obs, data$obs[[class_column]])
    
    # Initialize storage for indices
    selected_indices <- c()
    
    # First, ensure at least 2 samples per class
    for (class_data in class_groups) {
      if (nrow(class_data) >= 2) {
        selected_indices <- c(selected_indices, sample(rownames(class_data), size = 2, replace = FALSE))
      } else {
        stop("Insufficient samples to guarantee at least 2 per class in proportion: ", proportions[i])
      }
    }
    
    # Determine remaining samples to pick
    remaining_size <- size - length(selected_indices)
    
    if (remaining_size > 0) {
      # Remove already selected indices
      remaining_data <- data[!(rownames(data) %in% selected_indices), ]
      
      # Randomly sample remaining indices
      additional_indices <- sample(rownames(remaining_data), size = remaining_size, replace = FALSE)
      selected_indices <- c(selected_indices, additional_indices)
    }
    
    # Save the sampled subset
    samples[[sample_name]] <- data[selected_indices, ]
  }
  
  return(samples)
}

# Function to perform a stratified train/test split
train_test_split <- function(adata, train_size = 0.8, seed, class_column) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Extract the class labels from the AnnData object (assumed to be in metadata)
  class_labels <- adata$obs[[class_column]]
  
  # Ensure class_labels is a factor
  class_labels <- as.factor(class_labels)
  
  # Perform a stratified split based on the class_column
  train_indices <- NULL
  test_indices <- NULL
  
  # Stratified sampling: For each class, sample train/test indices proportionally
  for (class in levels(class_labels)) {
    # Get the indices of the samples belonging to the current class
    class_indices <- which(class_labels == class)
    
    # Number of samples in this class
    class_size <- length(class_indices)
    
    # Calculate the train and test sizes
    train_class_size <- floor(train_size * class_size)
    test_class_size <- class_size - train_class_size
    
    # Shuffle the class indices and split them
    set.seed(seed)  # Ensure reproducibility within each class
    shuffled_indices <- sample(class_indices, class_size, replace = FALSE)
    
    # Assign train/test indices
    train_indices <- c(train_indices, shuffled_indices[1:train_class_size])
    test_indices <- c(test_indices, shuffled_indices[(train_class_size + 1):(train_class_size + test_class_size)])
  }
  
  # Create the train and test AnnData objects
  train_adata <- adata[train_indices, ]
  test_adata <- adata[test_indices, ]
  
  # Return the train and test AnnData objects
  return(list(train = train_adata, test = test_adata))
}













# Functional analysis

# Function to read geneset libraries
read_enrichR_database <- function(path){
  # Read the file:
  res <- readLines(path)
  # Split each line by tabs:
  res <- strsplit(res, "\t")
  # First entry per tab is the gene set name:
  names(res) <- sapply(res, function(x) x[1])
  # Entries 3 to end are the gene names
  lapply(res, function(x) x[3:length(x)])
}

# Function to perform GSEA with fgsea on pvals and logFC
perform_gsea <- function(dge_results, database, rank_by = "logFC") {
  
  # Validate the rank_by argument
  if (!(rank_by %in% colnames(dge_results))) {
    stop("The specified rank_by column does not exist in the dge_results dataframe.")
  }
  
  # Rank genes by the specified rank_by column using the new pipe
  ranked_stats <- dge_results |>
    arrange(desc(.data[[rank_by]])) |>
    distinct(Feature, .keep_all = TRUE) |>
    pull(.data[[rank_by]], Feature)
  
  # Check if all values in ranked_stats are positive
  score_type <- if (all(ranked_stats > 0)) {
    "pos"  # Use "pos" if all values are positive
  } else {
    "std"  # Use "std" if there are negative values
  }
  
  # Run fgsea with the dynamic scoreType
  fgsea_results <- fgsea(pathways = database, stats = ranked_stats, scoreType = score_type) |>
    mutate(ranked_by = rank_by)
  
  return(fgsea_results)
}

# Quality control
classify_meta <- function(df, covariates) {
  # Initialize lists to store continuous, discrete, and boolean covariates
  classified_covariates <- list(
    continuous = c(),
    discrete = c(),
    boolean = c()
  )
  
  # Iterate through each covariate
  for (covariate in covariates) {
    # Skip covariates not present in the DataFrame
    if (!covariate %in% colnames(df)) {
      next
    }
    
    # Extract the column data
    col_data <- df[[covariate]]
    unique_values <- length(unique(col_data))
    
    # Check if column contains string representations of boolean values (TRUE/FALSE)
    if (is.character(col_data) || is.factor(col_data)) {
      # Check if all values are variations of "TRUE" or "FALSE"
      possible_boolean_values <- c("TRUE", "FALSE", "True", "False", "true", "false")
      if (all(col_data %in% possible_boolean_values)) {
        # Convert to logical type (TRUE/FALSE)
        col_data <- col_data == "TRUE" | col_data == "True" | col_data == "true"
        classified_covariates$boolean <- c(classified_covariates$boolean, covariate)
      } else {
        classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
      }
    } else if (is.logical(col_data)) {
      classified_covariates$boolean <- c(classified_covariates$boolean, covariate)
    } else if (is.factor(col_data)) {
      classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
    } else if (is.numeric(col_data)) {
      # Heuristic: Numerical columns with <= 10 unique values are discrete
      if (unique_values <= 10) {
        classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
      } else {
        classified_covariates$continuous <- c(classified_covariates$continuous, covariate)
      }
    } else {
      # Treat other data types as discrete by default
      classified_covariates$discrete <- c(classified_covariates$discrete, covariate)
    }
  }
  
  return(classified_covariates)
}


# PCA
# Function to plot PCA combinations for a single condition
plot_pca_for_condition <- function(pca_df, pca_result, condition, default_condition = NULL) {
  # Check if the condition column exists in pca_df
  if (!(condition %in% colnames(pca_df))) {
    stop(paste("Condition", condition, "not found in pca_df."))
  }
  
  # Dynamically detect the number of PCs in pca_df (assumes columns are named "PC1", "PC2", etc.)
  pc_columns <- grep("^PC[0-9]+$", colnames(pca_df), value = TRUE)
  if (length(pc_columns) < 2) {
    stop("At least two PCs are required for pairwise combinations.")
  }
  
  # Generate all combinations of the detected PCs
  pc_combinations <- combn(pc_columns, 2, simplify = FALSE)

  # Extract variance explained by each PC
  pc_variances <- summary(pca_result)$importance[2, ]
  
  # Create plots for each PC combination
  plot_list <- list()
  for (comb in pc_combinations) {
    x_pc <- comb[1]
    y_pc <- comb[2]

    # Calculate the variance explained by each PC in the combination
    x_var <- round(pc_variances[x_pc] * 100, 2)
    y_var <- round(pc_variances[y_pc] * 100, 2)
    
    # Dynamically handle the shape aesthetic if default_condition is not provided
    p <- ggplot(pca_df, aes(x = !!sym(x_pc), y = !!sym(y_pc), color = !!sym(condition))) +
      geom_point() +
      labs(
        x = paste0(x_pc, " (", x_var, "% Variance)"),
        y = paste0(y_pc, " (", y_var, "% Variance)")
      )
    
    # Add the shape aesthetic only if default_condition is provided
    if (!is.null(default_condition)) {
      if (!(default_condition %in% colnames(pca_df))) {
        stop(paste("Default condition", default_condition, "not found in pca_df."))
      }
      p <- p + aes(shape = !!sym(default_condition))
    }
  
    # Add the plot to the list
    plot_list[[paste(x_pc, y_pc, sep = "_")]] <- p
  }
  
  return(plot_list)
}


# Function to generate a list of t-SNE plots with different perplexities
generate_tsne_plots <- function(data, meta, condition, perplexity_range, default_condition = NULL, seed = 42) {
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Create an empty list to store plots
  plot_list <- list()
  
  # Loop through each value in the perplexity range
  for (perplexity in perplexity_range) {
    
    # Perform t-SNE with the current perplexity
    tsne_result <- Rtsne(data, dims = 2, perplexity = perplexity)
    
    # Create a data frame with t-SNE results and the condition labels
    tsne_data <- data.frame(
      TSNE1 = tsne_result$Y[, 1],
      TSNE2 = tsne_result$Y[, 2],
      Perplexity = factor(paste("Perplexity =", perplexity))
    )

    # Add meta
    tsne_data <- as_tibble(bind_cols(tsne_data, meta))
    
    # Create a ggplot for this perplexity value
    tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, 
                                       y = TSNE2, 
                                       color = !!sym(condition))) +
      geom_point() +
      facet_grid(cols = vars(Perplexity)) +
      labs(
        x = "t-SNE1",
        y = "t-SNE2"
      )
    
    # Add the shape aesthetic only if default_condition is provided
    if (!is.null(default_condition)) {
      if (!(default_condition %in% colnames(tsne_data))) {
        stop(paste("Default condition", default_condition, "not found in meta data."))
      }
      tsne_plot <- tsne_plot + aes(shape = !!sym(default_condition))
    }
  
    # Add the plot to the list
    plot_list[[paste("Perplexity", perplexity)]] <- tsne_plot
  }
  
  return(plot_list)
}
