# Standard libraries
import yaml
import pandas as pd
import numpy as np
import anndata
import pickle

# Tools for command line input
import sys
import argparse

# Import library
import importlib

# Subprocess to run R
import subprocess

# Parellel library
from joblib import Parallel, delayed

# Machine learning libraries
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

# Visualization
import seaborn as sns

# Custom functions
from py_utils import *

# Command line parser
parser = argparse.ArgumentParser()
parser.add_argument('settings_yml', type=str, help='path to setup file')

# Check if script is run from command line
if len(sys.argv) > 1:
    args = parser.parse_args()
    YAML_FILE = args.settings_yml
else:
    # Dev settings 
    YAML_FILE = 'dev_settings.yml'
    print('Running interactive mode for development.')

# Here starts the script   

# Initalize settings class
setup = Setup(YAML_FILE)
# Create output directory
setup.make_output_directory()
# Load settings
with open(setup.out('Settings.yml'), 'w') as f:
    yaml.dump(setup.final_config, f)
# Log
setup.log('Settings done')


# Error handler
error_file = setup.out("Error.log")
error_handler = ErrorHandler(log_file=error_file)

# Data loading
setup.log("Loading data")
data = anndata.read_h5ad(setup.data_path)


setup.log("Check data")
data_validator = DataValidator(data=data, setup=setup, error_handler=error_handler)
# Check if defined settings are present in the data
data_validator.data_settings_compatibility()
# Check for NAs in molecular data
data_validator.validate_na_counts()
# Check for negative counts in molecular data
data_validator.validate_negative_counts()



# Get European data
setup.log('Define ancestry')
eur_data = data[data.obs[setup.classification['ancestry_column']] == setup.classification['train_ancestry']]

# Define classification task 
setup.log('Define classification')
eur_data = (eur_data[eur_data.obs[setup.classification['output_column']]
                     .isin(setup.classification['comparison'])]
                     )

setup.log('Setting seed')
# Seed is used to generate different European subsets
seed = setup.seed
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# Randomly split EUR into two subsets (without stratification)
# 1. Sample random index to split the data frame (50 50)
# 2. Subset data based on these indexes 
split_index = eur_data.obs.sample(frac=0.5, random_state=seed).index
# a) The 'constant_split' is keept as ground_truth (no further splitting)
# b) The 'sampling_split' is used to randomly sample by size
#   i) The only prerequesite is that for each class there are at least two replicates
#      E.g. Cancer_1 = 2 observations, Cancer_2 = 2 observations
#      This is needed by limma eBayes
constant_split = eur_data[split_index, ]
sampling_split = eur_data[~eur_data.obs.index.isin(split_index), ]

# Save the index (observations) of the 'constant_split'
# Remaining Europeans
save_str = f'Obs_constant.yml'
with open(setup.out(save_str), 'w') as f:
    yaml.dump(constant_split.obs_names.to_list(), f)

# Further downsampling
# 1. Subset the 'sampling_split' by proportion
#   Prerequisite are at least 2 samples per class in each further subset
# 2. Get the observation names for each split
proportions = [1.0, 0.8, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02]
subset_list = sample_by_size(adata=sampling_split, 
                             props=proportions, 
                             output_column=setup.classification["output_column"],
                             seed=seed)

# Save all the observations to load in the Rscript
for subset in subset_list:
    # Save observation with the name of the subset
    save_str = f"Obs_{subset}.yml"
    observation = subset_list[subset].obs_names.to_list()
    with open(setup.out(save_str), 'w') as f:
        yaml.dump(observation, f)

# RUN R 
# ScriptRunner defines bin for R and python
R = ScriptRunner(r_path='/usr/bin/Rscript', py_path='/opt/conda/envs/ancestry/bin/python3')
settings_file = os.path.join(os.getcwd(), setup.output_directory, 'Settings.yml')

# Statistical analysis (DGE)
# Runs R script
# Uses files from the disk and writes to the disk
# All files are strored in specified 'output_directory'
setup.log('DGE')
R.run_script(script_path=os.path.join(os.getcwd(), 'eur_subsetting_dge.R'),
             args=[settings_file])

# Subset features that are used in DGE
# 1. Load Features
# 2. Subset 'constant_split'
# 3. Subset other subsets in 'subset_list'
p = os.path.join(setup.output_directory, 'Features.yml')
dge_features = yaml.safe_load(Path(p).read_text())

# 'constant_split'
constant_split = constant_split[:, dge_features]
# Other subsets
filtered_subset = []
for subset in subset_list:
    subset = subset_list[subset]
    subset = subset[:, dge_features]
    # Append
    filtered_subset.append(subset)

# Machine learning
setup.log('Sklearn data')
# Encoding for the 'output_column'
encoder = {}
for i,j in enumerate(setup.classification['comparison']):
    encoder.update({j: i})
inv_encoder = {v: k for k, v in encoder.items()}

# Transform anndata to np.array
# 'constant_split'
# 1. Transform data type
# 2. Normalize
constant_X = np.array(constant_split.X)
constant_X = normalize(eval(setup.normalization), constant_X)

constant_y = encode_y(constant_split.obs[setup.classification['output_column']], encoder)

# Other subsets
# 1. Iterate over all subsets
# 2. Transform data type
# 3. Normalize X
# 4. Append feature values to X and labels to y
subset_X = []
subset_y = []
for subset in filtered_subset:
    X = np.array(subset.X)
    # Normalization
    X = normalize(eval(setup.normalization), X)
    y = encode_y(subset.obs[setup.classification['output_column']], encoder)
    # Append
    subset_X.append(X)
    subset_y.append(y)

# Assertion:
# Check for dimensions
for i in zip(subset_X, subset_y):
   assert i[0].shape[0] == i[1].shape[0], \
    "Shapes of vectors do not align!"


print("Starting machine learning.")
setup.log("Machine learning")
# How many cpus per algorithm
total_jobs = setup.njobs  # Total jobs allowed
num_algorithms = len(setup.algorithms)
jobs_per_algorithm = total_jobs // num_algorithms
# Settings in pickable form (Setup class is non pickable)
pickable_settings = setup.return_settings()

# Training (parallel)
best_models = Parallel(n_jobs=num_algorithms)(
    delayed(train_algorithm)(algo_name, pickable_settings, 
                             constant_X, constant_y, 
                             jobs_per_algorithm)
    for algo_name in setup.algorithms
)

# Evaluation
metric_algo_list = []
for algo_name, best_model in zip(setup.algorithms, best_models):
    if best_model:
        metric_combined = evaluate_algorithm_eur_subsetting(algo_name, setup, best_model, subset_X, subset_y)
        metric_algo_list.append(metric_combined)
    else:
        print(f"Skipping evaluation for {algo_name} because training failed.")

# Join metric dataframes
merged_df = metric_algo_list[0]
# Iterate over the remaining dataframes and merge them
for df in metric_algo_list[1:]:
    merged_df = pd.merge(merged_df, df, on=["n_test_ancestry", "Metric"], how="outer")

# Add additional information
merged_df["Ancestry"] = setup.classification['train_ancestry'].upper()
merged_df["Seed"] = setup.seed
merged_df["n_train_ancestry"] = constant_split.shape[0]

# Save
merged_df.to_csv(setup.out('Metric_ml.csv'), index=False)



# # Iterate over different algortihms
# metric_algo_list = list()
# for algo_name in setup.algorithms:

#     setup.log(algo_name)
#     # Load the sklearn class
#     module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
#     algo_class = getattr(importlib.import_module(module_name), algo_name)
#     # Create instance of the algortihm
#     algo_instance = algo_class(random_state=setup.seed)

#     # Cross-validation for hyperparameter search
#     cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
#     # Load the grid_search parameters
#     search_space = setup.grid_search[algo_name]

#     setup.log("Training")
#     # Train with grid search and extract best model
#     grid_search = GridSearchCV(estimator=algo_instance, param_grid=search_space,
#                                cv=cv_hyp, scoring='f1_weighted', 
#                                refit=True, n_jobs=setup.njobs
#                                )
#     grid_search.fit(constant_X, constant_y)
#     best_m = grid_search.best_estimator_ 

#     # Save hyperparameters 
#     hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
#     hyp.to_csv(setup.out(f"Hyperparameters_{algo_name}.csv"), index=False)

#     # Save model
#     if setup.save_model:
#         with open(setup.out(f"{algo_name}_{setup.id}.pkl"), 'wb') as f:
#             pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

#     # Training finished
#     print(f"{algo_name} training done.")

#     # Test the model on all subsets
#     setup.log('Testing')
#     # Iterate over all X and y
#     # 1. Calculate propabilities
#     # 2. Calculate metric
#     # Initialize lists to store intermediate results
#     propability_list = []
#     metric_list = []
#     for X, y in zip(subset_X, subset_y):
#         # Get predicted probabilities
#         test_y_hat = pd.DataFrame(best_m.predict_proba(X))
#         test_y_hat["y"] = y
#         test_y_hat["n_test_ancestry"] = y.shape[0]
#         # Collect probabilities
#         propability_list.append(test_y_hat)

#         # Calculate metrics
#         test_probabilities, test_y = test_y_hat.iloc[:, 1].values, test_y_hat["y"].values
#         auc = roc_auc_score(y_true=test_y, y_score=test_probabilities)
#         metric = {algo_name: auc, "n_test_ancestry": y.shape[0], "Metric": "ROC_AUC"}
#         # Collect metrics
#         metric_list.append(metric)
    
#     # Validation finished
#     print(f"{algo_name} validation done")

#     # Combine probabilities and save 
#     propabilities = pd.concat(propability_list, ignore_index=True)
#     propabilities.to_csv(setup.out(f"Probabilities_{algo_name}.csv"), index=False)

#     # Combine metric and add algo information
#     metric_combined = pd.DataFrame(metric_list)

#     # Add to algo_list
#     metric_algo_list.append(metric_combined)

# # Join metric dataframes
# merged_df = metric_algo_list[0]
# # Iterate over the remaining dataframes and merge them
# for df in metric_algo_list[1:]:
#     merged_df = pd.merge(merged_df, df, on=["n_test_ancestry", "Metric"], how="outer")

# # Add additional information
# merged_df["Ancestry"] = setup.classification['train_ancestry'].upper()
# merged_df["Seed"] = setup.seed
# merged_df["n_train_ancestry"] = constant_split.shape[0]

# # Save
# merged_df.to_csv(setup.out('Metric_ml.csv'), index=False)
