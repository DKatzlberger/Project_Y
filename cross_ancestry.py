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
from sklearn.ensemble import RandomForestClassifier
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
    YAML_FILE = "dev_settings.yml"
    print("Running interactive mode for development.")

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


# Define training and inferred ancestry
setup.log('Define ancestries')
eur_data = data[data.obs[setup.classification['ancestry_column']] == setup.classification['train_ancestry']]
inf_data = data[data.obs[setup.classification['ancestry_column']] == setup.classification['infer_ancestry']]

# Define classification task 
setup.log('Define classification')
# Check that classification
eur_data = (eur_data[eur_data.obs[setup.classification['output_column']]
                     .isin(setup.classification['comparison'])]
                     )
inf_data = (inf_data[inf_data.obs[setup.classification['output_column']]
                     .isin(setup.classification['comparison'])]
                     )

# Assertion: Check min sample sizes per class
data_validator.check_min_samples_per_class(
    data=eur_data.obs,
    column=setup.classification['output_column'],
    min_samples=setup.sample_cutoff,
    data_name="European data"
    )

setup.log('Setting seed')
# Seed is used to generate different European subsets
seed = setup.seed
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)


setup.log("Create subset")
# Classify covariates
covariate_list = setup.classification.get("covariate", None)
if covariate_list:
    covariate_types = classify_covariates(eur_data.obs, covariates=covariate_list)
else:
    covariate_types = {'continuous': [], 'discrete': []}

# Stratify based on output column and discrete covariates
to_stratify_columns = [setup.classification["output_column"]] + covariate_types["discrete"]
# Create dictionary with frequencies
frequencies = inf_data.obs.groupby(to_stratify_columns, observed=False).size().to_dict()

# Split the training data into train and test set (Based on frequency of compared ancestry)
train_data, test_data = stratified_subset(eur_data, group_columns=to_stratify_columns, freq_dict=frequencies, seed=seed)
train_idx, test_idx, inf_idx = train_data.obs_names, test_data.obs_names, inf_data.obs_names

# Assertion: Check for data leakage
data_validator.check_data_leakage(train_idx, test_idx)
# Assertion: Check min sample size per class
# Train data
data_validator.check_min_samples_per_class(
    data=train_data.obs,
    column=setup.classification['output_column'],
    min_samples=2,
    data_name="Train data"
    )
# Test data
data_validator.check_min_samples_per_class(
    data=test_data.obs,
    column=setup.classification['output_column'],
    min_samples=2,
    data_name="Test data"
    )
# Inf data
data_validator.check_min_samples_per_class(
    data=inf_data.obs,
    column=setup.classification['output_column'],
    min_samples=2,
    data_name="Inf data"
    )

# TODO - Assertion: Covariate
# data_validator.check_covariate(
#     data=inf_data.obs,
#     column=setup.classification['output_column'],
#     min_samples=2,
#     data_name="Inf data"
#     )

# Save the samples that are used in train and test (use the same samples for 'stats')
# Remaining Europeans
save_str = f'Obs_train.yml'
with open(setup.out(save_str), 'w') as f:
    yaml.dump(train_idx.to_list(), f)

# European subset
save_str = f'Obs_test.yml'
with open(setup.out(save_str), 'w') as f:
    yaml.dump(test_idx.to_list(), f)

# Inferred Ancestry
save_str = f'Obs_inf.yml'
with open(setup.out(save_str), 'w') as f:
    yaml.dump(inf_idx.to_list(), f)

# RUN R define R to run
R = ScriptRunner(r_path='/usr/bin/Rscript', py_path='/opt/conda/envs/ancestry/bin/python3')
settings_file = os.path.join(os.getcwd(), setup.output_directory, 'Settings.yml')

# Statistical analysis (DGE)
setup.log('DGE')
script = "cross_ancestry_dge.R"
R.run_script(script_path=os.path.join(os.getcwd(), script),
             args=[settings_file])

# Subset features that are used in DGE
p = os.path.join(setup.output_directory, 'Features.yml')
dge_features = yaml.safe_load(Path(p).read_text())

train_data = train_data[:, dge_features]
test_data = test_data[:, dge_features]
inf_data = inf_data[:, dge_features]

# Machine learning
setup.log('Sklearn data')

# Encoding for the output column
encoder = {}
for i,j in enumerate(setup.classification['comparison']):
    encoder.update({j: i})
inv_encoder = {v: k for k, v in encoder.items()}

# Transform anndata to np.array
train_X = np.array(train_data.X)
train_y = encode_y(train_data.obs[setup.classification['output_column']], encoder)

test_X = np.array(test_data.X)
test_y = encode_y(test_data.obs[setup.classification['output_column']], encoder)

inf_X = np.array(inf_data.X)
inf_y = encode_y(inf_data.obs[setup.classification['output_column']], encoder)


# Normalization of features 
setup.log('Normalization')
train_X = normalize(eval(setup.normalization), train_X)
test_X = normalize(eval(setup.normalization), test_X) 
inf_X = normalize(eval(setup.normalization), inf_X) 

# Assertion: Check arrays
# Labels
assert_str = f'There have to be at least 2 classes, {np.unique(train_y).shape[0]} are given.'
assert np.unique(train_y).shape[0] >= 2, assert_str
# Dimension
assert_str = f'Dimensions of label vector ({train_y.shape[0]}) and feature matrix ({train_X.shape[0]}) do not align.'
assert train_y.shape[0] == train_X.shape[0], assert_str


print('Starting machine learning.')
setup.log('Machine learning')
# How many cpus per algorithm
total_jobs = setup.njobs  # Total jobs allowed
num_algorithms = len(setup.algorithms)
jobs_per_algorithm = total_jobs // num_algorithms
# Settings in pickable form (Setup class is non pickable)
pickable_settings = setup.return_settings()

# Training (parallel)
best_models = Parallel(n_jobs=num_algorithms)(
    delayed(train_algorithm)(algo_name, pickable_settings, 
                             train_X, train_y, 
                             jobs_per_algorithm)
    for algo_name in setup.algorithms
)

# Evaluation
y_hat_list = []
# Feature names for interpretations
feature_names = train_data.var_names
for algo_name, best_model in zip(setup.algorithms, best_models):
    if best_model:
          y_hat = evaluate_algorithm_cross_ancestry(
              algo_name, 
              best_model,
              setup, 
              test_X, 
              test_y, 
              inf_X, 
              inf_y, 
              feature_names,
              encoder = inv_encoder
              )
          # Append to list of probabilities
          y_hat_list.append(y_hat)
    else:
        print(f"Skipping evaluation for {algo_name} because training failed.")

# Combine probabilities across algorithms
y_hat = pd.concat(y_hat_list)
# Add information
y_hat["n_train_ancestry"] = len(train_idx)
y_hat["n_inf_ancestry"] = len(inf_idx)
# Save
y_hat.to_csv(setup.out(f"Probabilities.csv"), index=False)

# Calculate metric
# Multiclass calculation
if setup.classification['multiclass']:
    print('Score metrics for multiclass needs to be developed.')

# Binary calculation
else:
    R.run_script(script_path=os.path.join(os.getcwd(), 'binary_metric_calculation.py'),
                    args=[settings_file])

# Visualization
if setup.visual_val:
    setup.log('ML visualization')

    # Run the script to visualize
    R.run_script(script_path=os.path.join(os.getcwd(), 'run_visualization.R'), 
                 args=[settings_file])




# # Loop to train 
# for algo_name in setup.algorithms:

#     setup.log(algo_name)
#     # Load the sklearn class
#     module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
#     algo_class = getattr(importlib.import_module(module_name), algo_name)
#     # Create instance of the algortihm
#     algo_instance = algo_class(random_state=setup.seed)


#     # Cross-validation for hyperparameter search
#     cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=setup.seed)
#     # Load the grid_search parameters
#     search_space = setup.grid_search[algo_name]

#     setup.log("Training")
#     # Train with grid search and extract best model
#     grid_search = GridSearchCV(estimator=algo_instance, param_grid=search_space,
#                                cv=cv_hyp, scoring='f1_weighted', 
#                                refit=True, n_jobs=setup.njobs
#                                )
#     grid_search.fit(train_X, train_y)
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

#     setup.log("Testing")
#     # Validate prediction performance on EUR-subset
#     test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
#     test_y_hat["y"] = test_y
#     # Save propabilities
#     test_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_test.csv"), index=False)

#     setup.log("Infering")
#     # Validate prediction performance on ancestry
#     inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
#     inf_y_hat["y"] = inf_y
#     # Save propabilities
#     inf_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_inf.csv"), index=False)

#     # Validation finished
#     print(f"{algo_name} validation done")

#     setup.log("Model interpretations")
#     # Extract feature importances
#     feature_names = train_data.var_names
#     feature_importance = extract_feature_importance(best_m, feature_names)
#     # Save feature importances
#     feature_importance.to_csv(setup.out(f"Feature_importance_{algo_name}.csv"), index=False)