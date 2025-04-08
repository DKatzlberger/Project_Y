# Project_Y
Genotype-phenotype generalizability across ancestries

The repository (currently) contains:
1. A running script "interactions.R"

Input to the script:
1. Settings in yml format.



Output:
1. A directory with results.

"interactions.R" will create a directory within your specified `output_directory`. 
The name of this directory is structured as followed
```r
analysis_name    <- "interactions" # fixed
class_0          <- "class_0"      # modified by user
class_1          <- "class_1"      # modified by user
output_directory <- "development"  # modified by user

dir_name              <- paste0(class_0, "_vs_", class_1, "_", analysis_name)
path_to_save_location <- file.path(output_directory, dir_name)
```