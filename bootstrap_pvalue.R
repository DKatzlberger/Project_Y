# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    # Parallelization
    library(parallel)
    # Statistics 
    library(coin)
    library(fgsea)
    # Visualization
    library(patchwork)
    library(ggrepel)
    }
)
# Source custom functions
source("r_utils.R")
source("figure_themes.R")


comparison      <- "Basal_vs_non-Basal"
analysis_suffix <- "cross_ancestry"
match_pattern   <- paste0(comparison, ".*cross_ancestry_\\d+$")

# Load dirs with bootstraps metric
in_vscratch                 <- file.path("data", "runs1")
in_vscratch_dirs            <- list.dirs(in_vscratch, full.names = TRUE, recursive = FALSE)
in_vscratch_dirs_matched    <- grep(match_pattern, in_vscratch_dirs, value = TRUE)

# Output directory
out_vscratch <- file.path("data", "bootstrap_pvalue")
if (!dir.exists(out_vscratch)) {
  dir.create(out_vscratch, recursive = TRUE)
}


# Current filter
admix_dirs <- in_vscratch_dirs_matched[grep("ADMIX.*1$", in_vscratch_dirs_matched)]

# Load setting and metric files
train_all <- data.frame()
test_all  <- data.frame()
inf_all   <- data.frame()
for (dir in in_vscratch_dirs_matched){

    # Data path
    train_results   <- file.path(dir, "Limma_contrast_train.csv")
    test_results    <- file.path(dir, "Limma_contrast_test.csv")
    inf_results     <- file.path(dir, "Limma_contrast_inf.csv")

    # Load data
    if (file.exists(train_results)) {
        train_results_data <- read.csv(train_results)
    } else {
        message(paste("Train results file not found in:", dir))
    }
    
    if (file.exists(test_results)) {
        test_results_data <- read.csv(test_results)
    } else {
        message(paste("Test results file not found in:", dir))
    }
    
    if (file.exists(inf_results)) {
        inf_results_data <- read.csv(inf_results)
    } else {
        message(paste("Inference results file not found in:", dir))
    }

    # Seeting path
    settings_file_path <- file.path(dir, "Settings.yml")

    # Load settings
    if (file.exists(settings_file_path)) {
        # Load the settings from the YAML file
        settings_data <- yaml.load_file(settings_file_path)
    } else{
        message(paste("No settings file in,", dir))
    }

    # Currently seed info
    seed_info <- settings_data$seed

    print(seed_info)

    # Add info
    train_results_data$Seed <- seed_info
    test_results_data$Seed  <- seed_info
    inf_results_data$Seed   <- seed_info

    # Append
    train_all <- bind_rows(train_results_data, train_all)
    test_all  <- bind_rows(test_results_data, test_all)
    inf_all   <- bind_rows(inf_results_data, inf_all)
}

# Correlation
XY <- train_all |>
  select(Feature, logFC_X = logFC, Seed) |>
  inner_join(test_all |> select(Feature, logFC_Y = logFC, Seed), by = c("Feature", "Seed"))

XZ <- train_all |>
  select(Feature, logFC_X = logFC, Seed) |>
  inner_join(inf_all |> select(Feature, logFC_Z = logFC, Seed), by = c("Feature", "Seed"))

r_XY <- XY |>
  group_by(Seed) |>
  summarize(
    r_pearson   = cor(logFC_X, logFC_Y, method = "pearson"),
    r_spearman  = cor(logFC_X, logFC_Y, method = "spearman"),
    Prediction  = "Subset"
  )

r_XZ <- XZ |>
  group_by(Seed) |>
  summarize(
    r_pearson   = cor(logFC_X, logFC_Z, method = "pearson"),
    r_spearman  = cor(logFC_X, logFC_Z, method = "spearman"),
    Prediction  = "Ancestry"
  )

correlation <- 
  bind_rows(r_XY, r_XZ)

# Pivot
correlation_long <- correlation |>
  pivot_longer(
    cols = c(r_pearson, r_spearman), 
    names_to = "Metric", 
    values_to = "Value"
    ) |>
    mutate(
      Metric = recode(
        Metric, 
        r_pearson = "Pearson", 
        r_spearman = "Spearman"
        )
    )


# Plot
p <- correlation_long |>
  ggplot(
    aes(
      x    = Value,
      fill = Prediction
    )
  ) +
  geom_histogram(
    bins     = 100, 
    alpha    = 0.7, 
    position = "identity"
  ) +
  facet_grid(
    rows    = vars(Metric), 
    scales  = "free_y"
  ) +
  labs(
    y = "Count"
  ) +
  theme_nature_fonts(base_size = 8) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.position     = c(0.05, 0.95),   
    legend.justification = c(0, 1),        
  )

# Save
save_name <- file.path(out_vscratch, "Density.pdf")
save_ggplot(p, save_name, width = 4, height = 3)

# Bar plot
p <- correlation_long |>
  group_by(Metric, Prediction) |>
  summarize(
    mean_value = mean(Value),
    sd_value = sd(Value),
    .groups = "drop"
  ) |>
  ggplot(
    aes(
      x = fct_rev(Prediction),
      y = mean_value
      )
  ) +
  geom_bar(
    stat = "identity", 
    width = 0.7
  ) +
  geom_text(
    aes(
      label = round(mean_value, 3),  
      y = mean_value * 0.5  
    ),
    size = 2,  
    color = "white",  
  ) +
  geom_errorbar(
    aes(
      ymin = mean_value - sd_value, 
      ymax = mean_value + sd_value
      ), 
    width = 0.2, 
    position = position_dodge(0.7)
  ) +
  facet_grid(
    rows = vars(Metric)
  ) +
  labs(
    x = "Prediction",
    y = "Y"
  ) +
  theme_nature_fonts(base_size = 8) +
  theme_white_background() +
  theme_white_strip()

# Save
save_name <- file.path(out_vscratch, "Bar.pdf")
save_ggplot(p, save_name, width = 4, height = 3)



step_sizes <- seq(10, 100, by = 10)

plot_data <- map_df(step_sizes, function(n) {

    correlation_long |> 
      group_by(Metric, Prediction) |>
      slice_head(n = n) |>
      summarise(
        mean_corr = mean(Value),
        sd_corr = sd(Value),
        .groups = "drop"
      ) |>
      mutate(
        n_seeds = n
      )

  }
)




# Permutation test
permutation_test <- function(df, value, group, n_perm = 1000) {
  # Observed difference
  T_obs <- df |>
    group_by({{ group }}) |>
    summarise(mean = mean({{ value }}), .groups = "drop") |>
    summarise(diff = diff(mean)) |>
    pull(diff)

  perm_diffs <- numeric(n_perm)
  for (i in 1:n_perm) {
    shuffled <- df |>
      mutate(.group_shuffled = sample(pull(df, {{ group }}))) |>
      group_by(.group_shuffled) |>
      summarise(mean = mean({{ value }}), .groups = "drop") |>
      summarise(diff = diff(mean)) |>
      pull(diff)

    perm_diffs[i] <- shuffled
  }

  # One sided
  n_extreme <- sum(perm_diffs <= T_obs)
  p_value <- n_extreme / n_perm

  return(list(
    T_obs = T_obs,
    n_perm = n_perm,
    n_extreme = n_extreme,
    p_value = p_value,
    perm_diffs = perm_diffs
  ))
}


# Simulate p-value at increasing numbers of seeds
r_pearson <- correlation |>
  select(-r_spearman, -Seed)

step_sizes <- seq(10, 100, by = 10)

results <- map_df(step_sizes, function(n) {
  sampled_df <- r_pearson |>
    group_by(Prediction) |>
    slice_head(n = n) |>
    ungroup()
  
  # Run the permutation test and extract the p-value and differences
  perm_results <- permutation_test(
    sampled_df, 
    value = r_pearson, 
    group = Prediction, 
    n_perm = 10000
    )
  
  # Return a tibble of the results for the current step size
  tibble(
    n_seeds = n,
    p_value = perm_results$p_value,
    n_extreme = perm_results$n_extreme,
    T_obs = perm_results$T_obs,
    perm_diffs = list(perm_results$perm_diffs) # Store perm_diffs as list
  )
})

plot_data <- results |>
  unnest(cols = perm_diffs) |>
  mutate(n_seeds = factor(n_seeds))

seed_labels <- plot_data |>
  distinct(n_seeds, p_value) |>
  mutate(label = paste0(n_seeds, ", p = ", p_value)) |>
  select(n_seeds, label) |>
  deframe()

p <- ggplot(plot_data, aes(x = perm_diffs)) +
  geom_histogram(
    bins = 50, 
    fill = "skyblue", 
    color = NA
  ) +
  geom_vline(
    aes(xintercept = T_obs), 
    color = "red"
  ) +
  facet_grid(
    rows = vars(n_seeds), 
    labeller = labeller(n_seeds = as_labeller(seed_labels)),
    scales = "free"
  ) +  # Facet by n_seeds and allow y-axis scaling per facet
  labs(
    title = "Permutation Test: Is Ancestry Correlation Smaller?",
    subtitle = "Red Line: Observed mean difference",
    x = "Mean Difference (Shuffled)",
    y = "Frequency"
  ) +
  theme_nature_fonts(base_size = 10) +
  theme_white_background() +
  theme_white_strip()


# Save
save_name <- file.path(out_vscratch, "Permutation_test.pdf")
save_ggplot(p, save_name, width = 15, height = 15)



# Bootstrapping
step_sizes <- seq(10, 100, by = 10)

plot_data <- map_df(step_sizes, function(n) {
  
  # Subset for the current step size
  df <- correlation_long |>
    filter(Metric == "Pearson") |>
    group_by(Prediction) |>
    slice_head(n = n) |>
    ungroup() |>
    mutate(n_seeds = n)
  
  # Calculate observed test statistic: mean of Ancestry
  T_obs <- df |>
    filter(Prediction == "Ancestry") |>
    summarise(T_obs = mean(Value), .groups = "drop") |>
    pull(T_obs)

  # Extract EUR-subset values
  T_eur <- df |>
    filter(Prediction == "Subset") |>
    pull(Value)
  
  # One-sided p-value: is Subset >= Ancestry?
  p_value <- mean(T_eur >= T_obs)
  
  # Add stats back into df
  df |>
    mutate(
      T_obs = T_obs,
      p_value = p_value
    )
})

p_value_labels <- plot_data |>
  distinct(n_seeds, p_value) |>
  mutate(label = paste0(n_seeds, ", p = ", round(p_value, 4))) |>
  select(n_seeds, label) |>
  deframe()


p <- plot_data |>
  ggplot(
    aes(
      x    = Value,
      fill = Prediction
    )
  ) +
  geom_histogram(
    bins     = 100, 
    alpha    = 0.7, 
    position = "identity"
  ) +
  geom_vline(
    aes(xintercept = T_obs),
    color = "blue"
  ) +
  facet_grid(
    rows    = vars(n_seeds), 
    labeller = labeller(n_seeds = as_labeller(p_value_labels)),
    scales  = "free_y"
  ) +
  labs(
    y = "Count"
  ) +
  theme_nature_fonts(base_size = 8) +
  theme_white_background() +
  theme_white_strip() +
  theme_small_legend() +
  theme(
    legend.position     = c(0.05, 0.95),   
    legend.justification = c(0, 1),        
  )

# Save
save_name <- file.path(out_vscratch, "Bootstrapping.pdf")
save_ggplot(p, save_name, width = 5, height = 10)




r_mean_ancestry <- correlation |>
  select(-r_spearman) |>
  filter(Prediction == "Ancestry") |>
  summarise(
    r_pearson = mean(r_pearson),
    Prediction = "Ancestry"
  ) |>
  pull(r_pearson)

r_subset <- correlation |>
  select(-r_spearman) |>
  filter(Prediction == "Subset") |>
  pull(r_pearson)

# One sided
p_value <- mean(r_subset >= r_mean_ancestry)


