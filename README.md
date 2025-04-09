# Project_Y
Genotype-phenotype generalizability across ancestries

The repository (currently) contains:
1. A running script "interactions.R"

The script analysis which genes have interactions with ancestry in genotype-phenotype relationships.

**Input:**

***1. Settings:***

The script accepts one commandline argument, which is a path to a settings.yaml file. 
In the file you need to specify all required settings. In the snipped below some of the settings are displayed (more settings are required).
An example settings file for "interactions.R" can be found here [example_settings_interactions.yaml]().
```yaml
# Classification
class_0: your_class_0                  # Healthy, Cancer_1
class_1: your_class_1                  # Disease, Cancer_2
output_column: column_with_class       # Disease_status, Cancer_type
train_ancestry: your_train_ancestry    # EUR
infer_ancestry: your_infer_ancestry    # AFR
ancestry_column: column_with_ancestry  # Genetic, Self_reported
# Input
data_path: your/data/path.h5ad         # The script only works with .h5ad files
# Output
output_directory: your/save/location   # The script will create a directory at this place
```

***2. Data:***

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

**Output:**

"interactions.R" will create a directory at the place you specified in `output_directory`. 
This directory will contain the following files:

```
output_directory
|_ Settings.yaml
```