# Remove start up messages
suppressPackageStartupMessages(
    {
    # Standard libraries
    library(yaml)
    library(tidyverse)
    library(data.table)
    library(anndata)
    # Python
    library(reticulate)
    # Specify reticulate env
    use_condaenv("/opt/conda/envs/ancestry/bin/python")
    # Visualize
    library(patchwork)
    }
)

# Load data
adata <- read_h5ad("data/inputs/PanCanAtlas_COADREAD_RSEM.h5ad")
meta  <- adata$obs

# Visualize
cancer_type_detailed_plot <- meta |>
    ggplot(aes(x = cancer_type_detailed)) +
    geom_bar() +
    facet_grid(rows = vars(toupper(genetic_ancestry)), scale = "free") +
    labs(
        x = "Cancer type detailed",
        y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

ggsave(filename = "BRCA_cancer_type_detailed.pdf", 
       plot = cancer_type_detailed_plot, 
       path = "data/plots", 
       width = 3, height = 6.5
       )

subtype_plot <- meta |>
    ggplot(aes(x = subtype)) +
    geom_bar() +
    facet_grid(rows = vars(toupper(genetic_ancestry)), scale = "free") +
    labs(
        x = "Subtype",
        y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

ggsave(filename = "BRCA_subtype.pdf", 
       plot = subtype_plot, 
       path = "data/plots", 
       width = 3, height = 6.5
       )

# Patchwork
combined_plot <- cancer_type_detailed_plot + subtype_plot +
    plot_annotation(title = "BRCA cancer subtypes", 
                    theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave(filename = "BRCA_possible_comparisons.pdf", 
       plot = combined_plot, 
       path = "data/plots", 
       width = 5, height = 7
       )

meta <- fread("data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv")
pan_can_atlas_studies <- meta |>
    ggplot(
        aes(
            x = studyId
        )
    ) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

ggsave(filename = "PanCanAtlas_studies.pdf", 
       plot = pan_can_atlas_studies, 
       path = "data/plots", 
       width = 5, height = 7
       )