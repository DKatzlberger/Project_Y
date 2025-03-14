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
from sklearn.metrics import roc_auc_score, accuracy_score

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
setup.log('Settings done')

# Seed
setup.log("Setting seed")
seed = setup.seed
np.random.seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)

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
setup.log("Define ancestry")
eur_data = data[data.obs[setup.classification["ancestry_column"]] == setup.classification['train_ancestry']]

# Define classification task 
setup.log("Define classification")
eur_data = (eur_data[eur_data.obs[setup.classification["output_column"]].isin(setup.classification["comparison"])])

setup.log("Feature selection")
# Filtering based on eur_data
if setup.filter_features and setup.data_type == "expression":
    # Filter by expression
    cpm_data = cpm(eur_data.X)
    min_counts = min_counts_by_percentile(cpm_data, percentile=25)
    filtered_features = filter_by_expression(cpm_data, min_counts = min_counts)

    # Subset features
    eur_data = eur_data[:, filtered_features]

# Save features (for DGE)
with open(setup.out("Features.yml"), "w") as f:
    yaml.dump(eur_data.var_names.to_list(), f)


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

# Assertion: Check min sample sizes per class
data_validator.check_min_samples_per_class(
    data=constant_split.obs,
    column=setup.classification["output_column"],
    min_samples=setup.sample_cutoff,
    data_name="All train data"
    )

# Save the index (observations) of the 'constant_split'
# Remaining Europeans
save_str = f"Obs_constant.yml"
with open(setup.out(save_str), 'w') as f:
    yaml.dump(constant_split.obs_names.to_list(), f)

# Further downsampling
# 1. Subset the 'sampling_split' by proportion
#   Prerequisite are at least 2 samples per class in each further subset
# 2. Get the observation names for each split
proportions = [1.0, 0.8, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02]
subset_list = sample_by_size(
    adata=sampling_split, 
    props=proportions, 
    output_column=setup.classification["output_column"],
    seed=seed
    )


# Save all the observations
for subset in subset_list:
    # Assertion: Check for enough samples
    data_validator.check_min_samples_per_class(
        data=subset_list[subset].obs,
        column=setup.classification["output_column"],
        min_samples=2,
        data_name=f"{subset}"
    )

    # Save
    save_str = f"Obs_{subset}.yml"
    observation = subset_list[subset].obs_names.to_list()
    with open(setup.out(save_str), 'w') as f:
        yaml.dump(observation, f)

# Script runner
R = ScriptRunner(r_path="/usr/bin/Rscript", py_path='/opt/conda/envs/ancestry/bin/python3')
settings_file = os.path.join(os.getcwd(), setup.output_directory, "Settings.yml")

# Statistical analysis (DGE)
setup.log("DGE")
R.run_script(script_path=os.path.join(os.getcwd(), "eur_subsetting_dge.R"), args=[settings_file])

# Machine learning
setup.log("Sklearn data")
# Encoding for the output column
comparison = setup.classification["comparison"]
encoder = {j: i for i, j in enumerate(comparison)}
inv_encoder = {v: k for k, v in encoder.items()}

# Transform anndata to np.array
output_column = setup.classification["output_column"]
constant_X, constant_y = np.array(constant_split.X), encode_y(constant_split.obs[output_column], encoder)
# Normalization
setup.log("Normalization")
# Select normalization method
data_type = setup.data_type
ml_normalization = setup.ml_normalization
normalization_method = ml_normalization_methods[data_type][ml_normalization]
# Normalization
constant_X = normalization_method(constant_X)

# Test subsets
subset_X = []
subset_y = []
for sub in subset_list:
    subset = subset_list[sub]
    X = np.array(subset.X)
    # Normalization
    X = normalization_method(X)
    y = encode_y(subset.obs[output_column], encoder)
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
    delayed(train_algorithm)(
        algo_name, 
        pickable_settings, 
        constant_X, 
        constant_y, 
        jobs_per_algorithm
        )
    for algo_name in setup.algorithms
)

# Evaluation
probabilities_list = []
for algo_name, best_model in zip(setup.algorithms, best_models):
    if best_model:
        y_hat = evaluate_algorithm_eur_subsetting(
            algo_name, 
            best_model, 
            subset_X, 
            subset_y
            )
        # Append across algorithms
        probabilities_list.append(y_hat)
    else:
        print(f"Skipping evaluation for {algo_name} because training failed.")

# Combine across algorithm
probabilities = pd.concat(probabilities_list)
# Add additional information
probabilities = probabilities.rename(columns=inv_encoder)
probabilities["Seed"] = setup.seed
probabilities["Ancestry"] = setup["classification"]["train_ancestry"].upper()
# Save
probabilities.to_csv(setup.out(f"Probabilities.csv"), index=False)

# Claculate metric
metric_across_algo = []
for (algo_name, subset), algo_probabilities in probabilities.groupby(["Algorithm", "n_test_ancestry"]):
    # Extract true labels and probabilities
    subset_y, subset_probabilities = algo_probabilities.iloc[:, 2].values, algo_probabilities.iloc[:, 1].values
    # ROC AUC
    auc_score = roc_auc_score(y_true=subset_y, y_score=subset_probabilities)
    # Accuracy
    threshold = 0.5
    subset_preds = (subset_probabilities >= threshold).astype(int)
    acc_score = accuracy_score(y_true=subset_y, y_pred=subset_preds)
    # Prepare data
    df_data = {
        "Algorithm": algo_name, 
        "n_test_ancestry": subset,
        "ROC_AUC": auc_score,
        "Accuracy": acc_score,
        "Status": "Test", 
        "Prediction": "Subset",
        "Threshold": threshold
        }
    # Metric dataframe
    metric_df = pd.DataFrame(df_data, index=[0])
    # Add to across algorithms
    metric_across_algo.append(metric_df)

# Join metric dataframes
metric_across_algo = pd.concat(metric_across_algo)
# Add additional information
metric_across_algo["Seed"] = setup.seed
metric_across_algo["Ancestry"] = setup["classification"]["train_ancestry"].upper()

# Save
output_directory = setup["output_directory"]
metric_across_algo.to_csv(os.path.join(output_directory, "Contrast_metric_ml.csv"), index=False)

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
