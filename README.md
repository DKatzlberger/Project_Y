# Genotype-phenotype generalizability across ancestries

The repository (currently) contains:
1. A running script "interactions.R"

The script analysis which genes have interactions with ancestry in genotype-phenotype relationships. It uses the edgeR pipeline to compare each gene individually. The design of the model first fits a means model for the specified groups and afterwards fits a contrast to compare these groups. The script does the following comparison:
1. class_0 train_ancestry vs infer_ancestry (Baseline 1)
2. class_1 train_ancestry vs infer_ancestry (Baseline 2)
3. class_0 vs class_1 train_ancestry (Relationship 1)
4. class_0 vs class_1 infer_ancestry (Relationship 1)
5. class_0 vs class_1 train_ancestry vs infer_ancestry (Interaction)

## Input

***1. Settings:***

The script accepts one commandline argument, which is a path to a settings.yaml file. 
In the file you need to specify all required settings. An example settings file for "interactions.R" can be found here [example_settings_interactions.yaml](https://github.com/DKatzlberger/Project_Y/blob/main/example_settings_interactions.yaml). The snipped below shows all required settings:
```yaml
# Required settings
class_0: your_class_0                  
class_1: your_class_1                  
output_column: column_with_class       
train_ancestry: your_train_ancestry    
infer_ancestry: your_infer_ancestry    
ancestry_column: column_with_ancestry  
data_path: your/data/path.h5ad         # The script only works with .h5ad files
tech: omics_type                       
normalization: method
output_directory: your/save/location   # The script will create a directory at this place
```
Currently the script supports these omics and normalization methods:
| Omics type (tech)    | Supported normalization methods (normalization)                      |
|----------------------|----------------------------------------------------------------------|
| transcriptomics      | `"limma_voom"` `"normalize_log"` `"normalize_zscore"` `"raw"`        |
| methylation          | `"beta_to_mvals"` `"normalize_log"` `"normalize_zscore"` `"raw"`     |
| proteomics           | `"raw"`                                                              |

The script will substitute not specified settings with default settings. Default settings are specified in [default_settings_interactions.yaml](https://github.com/DKatzlberger/Project_Y/blob/main/default_settings_interactions.yaml) and can be modified.

***2. Data:***

The script handles count matrix with one value per gene.
The data is structured in an [AnnData](https://anndata.readthedocs.io/en/stable/) object. 
*X* is required to contain your molecular data in format *observations* x *features*.
*obs* (meta data) is required to include the `output_column` with specified `class_0`, `class_1` and it requires to contain the `ancestry_column` with specified `train_ancestry`, `infer_ancestry`.
```r
data_path <- "your/data/path.h5ad"
adata     <- read_h5ad(data_path)
adata
```
```
AnnData object with n_obs × n_vars = 794 × 20338
    obs: column_with_class, column_with_ancestry
```

## Output

"interactions.R" will create a directory at the place you specified in `output_directory`. 
This directory will contain the following files:

```
output_directory/
├── Features.yaml                      # Analyzed features
├── Settings.yaml                      # Used settings
├── Limma_means.csv                    # Result means model
├── Limma_contrast.csv                 # Results contrast
├── QC_mean_variance_trend.pdf         
├── QC_desnity_normalized_values.pdf 
├── QC_qq_normalized_values.pdf 
├── QC_sample_sizes.pdf
└── Visual_val/
    ├── Baseline_ma.pdf
    ├── Baseline_volcano.pdf
    ├── Relationship_ma.pdf
    ├── Relationship_volcano.pdf
    ├── Interaction_ma.pdf
    ├── Interaction_volcano.pdf
    ├── Interaction_heatmap.pdf
    └── Interaction_boxplot.pdf
```

## How to run the script
To run the analysis using a Singularity container, use the following bash script. 
This script executes the R script with a specified YAML configuration file inside the container:
```bash
singularity exec <SIF_PATH> Rscript <SCRIPT_PATH< <YAML_CONFIG>
```