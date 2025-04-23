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
source("figure_themes.R")

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


# Load samples
samples   <- read.table("data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients/GSE225846_series_matrix.txt", sep = "\t", comment.char = "!", header = TRUE)[-1]
n_samples <- length(samples)
# Parse and split lines
lines <- readLines("data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients/GSE225846_series_matrix.txt")
lines <- strsplit(lines, "\t")
lines <- lines[sapply(lines, length) > n_samples]
# Clean lines
clean_lines    <- lapply(lines, function(x) gsub('^\"|\"$', '', x))
filtered_lines <- clean_lines[sapply(clean_lines, function(x) grepl("^!", x[1]))]
# Named list
column_data <- lapply(filtered_lines, function(x) {
  column <- sub("^!", "", x[1])  
  values <- x[-1]  
  setNames(list(values), column)  
})
# Flatten 
flat_data <- unlist(column_data, recursive = FALSE)
data      <- as.data.frame(flat_data)
colnames(data)

# Important information
char_data <- data[, grep("^Sample_characteristics_ch1", colnames(data))]
new_cols  <- sapply(char_data[1, ], function(x) sub(":.*", "", x))
# Rename colnames
colnames(char_data) <- new_cols
colnames(char_data) <- gsub(" ", "_", colnames(char_data))
# Remove from value
char_data[]         <- lapply(char_data, function(col) sub("^[^:]*: ?", "", col))
# Add additional inforamtion
char_data["sample_title"]  <- gsub("^S_", "", data$Sample_title)
char_data["geo_accession"] <- data["Sample_geo_accession"]
char_data["tech"]          <- data["Sample_library_source"] 
char_data["strategy"]      <- data["Sample_library_strategy"] 
# Save
path_to_save_location <- "data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients"
save_name             <- file.path(path_to_save_location, "Meta_data.csv")
fwrite(char_data, save_name)

# Plot
p <- char_data |>
  ggplot(
    aes(
      x    = type,
      fill = race
    )
  ) +
  geom_bar(
    # position = "fill"
  ) +
  coord_flip() +
  labs(
    title = "Breast tumor and non-cancerous adjacent tissue ",
    x     = "Disease status",
    y     = "Count",
    fill  = "Self-reported ancesetry"
  )
# Save 
save_name <- file.path(path_to_save_location, "Meta_plot.pdf")
save_ggplot(p, save_name, width = 6, height = 3)

# Molecular data
data <- read_tsv("data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients/GSE225846_RawCountFile_rsemgenes.txt", col_names = TRUE)
data <- separate(data, gene_id, into = c("entrez", "symbol"), sep = "_", extra = "merge")
data <- select(data, -entrez)
# Unique gene id
unique_genes <- distinct(data, symbol, .keep_all = TRUE)
unique_genes <- unique_genes |>
  column_to_rownames("symbol") |>
  t() |>
  as.data.frame() |>
  rownames_to_column("sample_title")
# Save
save_name <- file.path(path_to_save_location, "Molecular_data.csv")
fwrite(unique_genes, save_name)

