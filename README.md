# Project_Y
Genotype-phenotype generalizability across ancestries

The repository (currently) contains:
1. A running script "interactions.R"

The script analysis which genes have interactions with ancestry in genotype-phenotype relationships.

Input:
1. Settings in yaml format
2. Data in h5ad format


1. Settings:

The script accepts one commandline argument, which is a `settings.yaml`file. 
In the file you specify all your settings:
```yaml
# Classification
class_0: your_class_0
class_1: your_class_1
output_column: column_with_class
train_ancestry: your_train_ancestry
infer_ancestry: your_infer_ancestry
ancestry_column: column_with_ancestry
# Input
data_path: your/data/path.h5ad          # The script only works with .h5ad files
```


Output:
1. A directory with results.

"interactions.R" will create a directory within your specified `output_directory`. 
The name of this directory is structured as followed
```r
class_0          <- "class_0"      # modified by user
class_1          <- "class_1"      # modified by user
train_ancestry   <- "ancestry_0"   # modified by user
infer_ancestry   <- "ancestry_1"   # modified by user
output_directory <- "development"  # modified by user

analysis_name         <- "interactions"
phenotypes            <- paste0(class_0, "_vs_", class_1)
ancestries            <- paste0(train_ancestry, "_to_", infer_ancestry)
dir_name              <- paste(c(phenotypes, ancestries, analysis_name), collapse = "_")
path_to_save_location <- file.path(output_directory, dir_name)

# > path_to_save_location
# [1] "development/class_0_vs_class_1_ancestry_0_to_ancestry_1_interactions"
```