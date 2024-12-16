suppressPackageStartupMessages({
    # Standard libraries
    library(tidyverse)
    library(data.table)
    # Visualization
    library(patchwork)
})

# Script to visulaize meta data downloaded from the cBioPortal
meta <- fread('data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv')

# Plot 'cancer_type_detailed' for each ancestry
meta_plot <- meta |> 
  as_tibble() |>
  filter(cancer_type_detailed != '') |>
  mutate(tumor_type = toupper(tumor_type),
         consensus_ancestry = toupper(pooled_consensus_ancestry)) |> 
  ggplot(aes(x = cancer_type_detailed)) +
  geom_bar(stat = 'count', aes(fill = after_stat(count) > 25)) +
  facet_grid(cols = vars(tumor_type),
             rows = vars(pooled_consensus_ancestry), 
             scales = 'free',
             space = 'free_x') +
  scale_fill_manual('Cancers colored by sample sizes',
                    values = c("#FC2D00","#008EFC"),
                    labels = c('TRUE'='> 25','FALSE'='< 25')) +
  xlab('Cancer type (detailed)') +
  ylab('Count') +
  ggtitle("TCGA's sub-cancer types") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        strip.text.x = element_text(size = 7),
        panel.spacing.x = unit(0.2, 'lines'),
        legend.position = 'bottom',
        legend.title.position = 'top',
        legend.title = element_text(hjust = 0.5)
        )

# Save the plot
ggsave('data/plots/Possible_prediction_tasks.pdf', width = 20, height = 15, scale = 1)
