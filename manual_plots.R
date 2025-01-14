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
    }
)

# Load data
adata <- read_h5ad("data/inputs/PanCanAtlas_BRCA_transcriptome_RSEM.h5ad")
meta  <- adata$obs
meta |>
    ggplot(aes(x = cancer_type_detailed)) +
    geom_bar() +
    facet_grid(rows = vars(genetic_ancestry), scale = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

unique(adata$obs$genetic_ancestry)
