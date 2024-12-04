# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Statistics 
    library(coin)
    }
)
# Source custom functions
source("r_utils.R")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the script is run interactively or with command-line arguments
if (length(args) > 0) {
  yaml_file <- args[1] 
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)
  # When the script is run from the command line then 'output_directory' is given
  # The pattern to extract all matchin directories is extracted from 'output_directory'
  output_path = setup$output_directory
  match_pattern <- gsub("_\\d+$", "", output_path)
  # TODO - Check if regular expression also works for other approaches

} else {
  print("Running interactive mode for development.")
  # Yaml file used for development (often an actual job script)
  yaml_file <- "job_settings.yml"
  # 1. Check if it's a valid yaml file
  # 2. Load the yaml file
  is_yml_file(yaml_file)
  setup <- yaml.load_file(yaml_file)

  # When run for development no 'output_directory' is specified
  # Hence, the 'output_directory' has to be contructed like in the 'process_yaml_... .py'

  # Construction:
  # Vscratch_dir is the place where the files are stored
  vscratch_dir_in = "data/runs"
  # Tag is used to specify which data it is e.g. TCGA, NIAGADS
  tag <- setup$tag
  # Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
  comparison <- paste0(tag, "_", paste(setup$classification$comparison, collapse = "_vs_"))
  # Analysis_name specifies what was analysed
  # E.g. comparsion of ancestries EUR_to_AMR, subsetting_EUR
  # This often is modified depending which analysis
  analysis_name  <-  "EUR_subsetting"
  # Combination of components to create the 'match_pattern'
  # The 'match_pattern' is used as pattern to extract all folders in the vscratch dir
  match_pattern <- paste0(comparison, "_", analysis_name)
}

# Extracting all folders in the 'vscratch_dir' that match 'match_pattern'
# 1. List all folders in 'vscratch_dir_in'
# 2. With 'match_pattern' extract matching folders
all_vscratch_dir_in <- list.dirs(vscratch_dir_in, full.names = TRUE, recursive = FALSE)
match_vsratch_dir <- grep(match_pattern, all_vscratch_dir_in, value = TRUE)

# Check if there were matching folders
if (length(match_vsratch_dir) == 0) {
  message("No matching folders found.")
}

# Save the results of the analysis: 'EUR_subsetting'
# 'Vscratch_dir_out' where summarized analysis are stored
vscratch_dir_out  <- "data/combined_runs"
path_to_save_location <- file.path(vscratch_dir_out, comparison, analysis_name)
# Create the directory also parents (if it does not exist)
# Create directory if it does not exists
if (!dir.exists(path_to_save_location)) {
  dir.create(path_to_save_location, recursive = TRUE)
}


# Loading and combining the matched folders
metric_all <- data.frame()
for (folder in match_vsratch_dir){
    # Generate path to csv_file
    # All csv_files are stored under 'Metric_subsetting.csv'
    # 1. Load the files
    # 2. Combine them into single data frame
    path_to_metric <- file.path(folder, "Metric_subsetting.csv")
    metric <- fread(path_to_metric)
    metric_all <- bind_rows(metric_all, metric) 
}

# Check for outliers
# TODO - Check for outliers per proportion of the subset using hampel filer 

# Summarize the across seeds
# Calculate mean, sd, sem 
summarized_metric <- metric_all |> 
    pivot_longer(
        cols = c(Pearson, Spearman),
        values_to = "Value",
        names_to = "Metric"
    ) |> 
    # Summarize by Seed 
    # Dont use Seed in 'group_by'
    group_by(Proportion_of_constant_split, n, Metric) |>  
    summarize(
        n_seeds = n(),
        mean_value = mean(Value, na.rm = TRUE),
        sd_value = sd(Value, na.rm = TRUE),
        se_value = sd(Value, na.rm = TRUE) / sqrt(n())
  ) |>
  mutate(
    n = as.factor(n),
    Proportion_of_constant_split = as.factor(Proportion_of_constant_split)
    )

# Plot the correlation
# Correlation axis from 0 to 1
common_y <- scale_y_continuous(
    limits = c(0, 1.1), 
    breaks = c(0, 0.5, 1))

# TODO 
test <- summarized_metric |>
  ggplot(
    aes(
        x = fct_rev(n),
        y = mean_value
    )
  ) +
  geom_bar(
        stat = "identity", 
        width = 0.7
  ) +
  geom_errorbar(
        aes(
            ymin = mean_value - sd_value, 
            ymax = mean_value + sd_value
            ), 
        width = 0.2, position = position_dodge(0.7)
  ) +
  common_y +
  facet_grid(
    rows = vars(Metric)
  ) +
  labs(
    # title = paste(gsub("_", " ", comparison)),
    x = "EUR subsets (n)",
    y = "Correlation coefficient"
  ) +
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(path_to_save_location, "Subsets.pdf"),
        plot = test, height = 6, width = 12)
