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

# Function to check if covariate is present in settings
get_covariate <- function(setup) {
  # Function to extract the covariate if present and valid
  if ("covariate" %in% names(setup$classification) && 
      !is.null(setup$classification$covariate) &&
      any(nzchar(setup$classification$covariate))) {
    return(setup$classification$covariate)
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

# Calculate Kolmogorow-Smirnow
ks_test <- function(data, group_col, value_col, group1 = "Test", group2 = "Inference") {
  
  # Filter the data for the two groups and pull the corresponding values
  group_data <- data |>
    filter(!!sym(group_col) %in% c(group1, group2)) |>
    select(!!sym(group_col), !!sym(value_col))
  
  # Extract the values for the two groups to compare
  group1_values <- group_data |> filter(!!sym(group_col) == group1) |> pull(!!sym(value_col))
  group2_values <- group_data |> filter(!!sym(group_col) == group2) |> pull(!!sym(value_col))
  
  # Perform Kolmogorov-Smirnov test
  ks_result <- ks.test(group1_values, group2_values)
  
  # Return the result of the KS test
  return(ks_result)
}

# Permutation test
permutation_test <- function(data, value_col, group_col) {
  
  # Check if any group has identical values (zero variance)
  var_check <- data |>
  group_by(!!sym(group_col)) |>
  summarise(var_value = var(!!sym(value_col)), .groups = 'drop')
  
  # If any group has zero variance, return p-value of 1
  if (any(var_check$var_value == 0)) {
    message("One or more groups have identical values. Returning p-value of 1.")
    return(NA)  # Return the p-value of NA directly
  }

  # Quote the column names for later use in evaluation
  value_col <- sym(value_col)  # Quoting the value column
  group_col <- sym(group_col)  # Quoting the group column
  
  # Prepare the data by selecting relevant columns and ensuring the group variable is a factor
  perm_data <- data |>
    select(!!value_col, !!group_col) |>
    mutate(!!group_col := as.factor(!!group_col))  # Ensure group_col is a factor
  
  # Construct the formula for the independence test
  formula <- as.formula(paste(deparse(value_col), "~", deparse(group_col)))
  
  # Perform the permutation test (independence test)
  perm_test_result <- independence_test(formula = formula, data = perm_data)
  
  # Return the result of the permutation test
  p_value = as.numeric(pvalue(perm_test_result))
  return(p_value)
}

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
perform_gsea <- function(dge_results, database_path, rank_by = "logFC") {
  
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
  
  # Read the database (path to file)
  db <- read_enrichR_database(database_path)
  
  # Run fgsea with the dynamic scoreType
  fgsea_results <- fgsea(pathways = db, stats = ranked_stats, scoreType = score_type) |>
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
plot_pca_for_condition <- function(pca_df, pca_result, condition, default_condition) {
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
    
    # Create a ggplot for the given PC combination
    p <- ggplot(pca_df, aes(x = !!sym(x_pc), y = !!sym(y_pc), 
                            color = !!sym(condition),
                            shape = !!sym(default_condition)
                            )
      ) +
      geom_point(size = 3) +
      labs(
        x = paste0(x_pc, " (", x_var, "% Variance)"),
        y = paste0(y_pc, " (", y_var, "% Variance)")
      )
  
    # Add the plot to the list
    plot_list[[paste(x_pc, y_pc, sep = "_")]] <- p
  }
  
  return(plot_list)
}


# Function to generate a list of t-SNE plots with different perplexities
generate_tsne_plots <- function(data, meta, condition, perplexity_range, default_condition, seed = 42) {
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
                                      color = !!sym(condition),
                                      shape = !!sym(default_condition))
      ) +
      geom_point() +
      facet_grid(cols = vars(Perplexity))
      labs(
        x = "t-SNE1",
        y = "t-SNE2"
      ) 
    
    # Add the plot to the list
    plot_list[[paste("Perplexity", perplexity)]] <- tsne_plot
  }
  
  return(plot_list)
}