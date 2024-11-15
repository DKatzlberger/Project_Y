library(tidyverse)
library(data.table)

# Path to the metric
directory_path <- "data/combined_runs"

# List all CSV files in the directory
csv_files <- list.files(
    path = directory_path, 
    pattern = "TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_.*_combined_ml_metric.csv", 
    full.names = TRUE)

# Read and combine data
combined_data <- csv_files |> 
    lapply(function(file) {
        df <- read.csv(file) 
        return(df)
    }) |> 
    bind_rows()

# Plot
combined_data |> 
    ggplot(
        aes(
            x = fct_rev(Prediction),
            y = ROC_AUC
        )
    ) +
    geom_col() +
    facet_grid(
        cols = vars(Ancestry)
    ) +
    xlab('Prediction') +
    ylab('ROC AUC') +
    ggtitle(combined_data$Comparison)

# Save the plot
base_dir <- 'data/plots'
comparison <- 'TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_ml_performance.pdf'
save_path <- file.path(base_dir, comparison)
ggsave(save_path)


# Differential gene expression
# List all CSV files in the directory
csv_files <- list.files(
    path = directory_path, 
    pattern = "TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_.*_combined_dge_metric.csv", 
    full.names = TRUE)

# Read and combine data
combined_data <- csv_files_1 |> 
    lapply(function(file) {
        df <- read.csv(file) 
        return(df)
    }) |> 
    bind_rows()

# Plot
combined_data |> 
    ggplot(
        aes(
            x = fct_rev(Prediction),
            y = Pearson_correlation
        )
    ) +
    geom_col() +
    facet_grid(
        cols = vars(Ancestry)
    ) +
    xlab('Prediction') +
    ylab('Pearson') +
    ggtitle(combined_data$Comparison)

# Save the plot
base_dir <- 'data/plots'
comparison <- 'TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_dge_performance.pdf'
save_path <- file.path(base_dir, comparison)
ggsave(save_path)
