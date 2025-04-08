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
from sklearn.metrics import roc_auc_score, accuracy_score

# Visualization
import seaborn as sns

# Custom functions
from py_utils import *

# Command line parser
parser = argparse.ArgumentParser()
parser.add_argument("settings_yml", type=str, help="path to setup file")

# Check if script is run from command line
if len(sys.argv) > 1:
    args = parser.parse_args()
    YAML_FILE = args.settings_yml
else:
    # Dev settings 
    # BRCA
    # data/inputs/settings/PanCanAtlas_BRCA_BETA_Basal_vs_non-Basal_EUR_to_ADMIX.yml

    # UCEC
    # data/inputs/settings/PanCanAtlas_UCEC_RSEM_CN-high_vs_non-CN-high_EUR_to_ADMIX.yml

    YAML_FILE = "dev_settings.yml"
    print("Running interactive mode for development.")

# Here starts the script   

# Initalize settings class
setup = Setup(YAML_FILE)
# Create output directory
setup.make_output_directory()

# Save settings
with open(setup.out("Settings.yml"), 'w') as f:
    yaml.dump(setup.final_config, f)
setup.log("Settings done")

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

# Validate data
setup.log("Check data")
data_validator = DataValidator(data=data, setup=setup, error_handler=error_handler)
# Check if defined settings are present in the data
data_validator.data_settings_compatibility()
# Check for NAs in molecular data
data_validator.validate_na_counts()
# Check for negative counts in molecular data
data_validator.validate_negative_counts()


# Define classification task 
setup.log("Define classification")
data = (data[data.obs[setup.output_column].isin(setup.comparison)])

setup.log("Feature selection")
# Filtering based on all of the data (consistent with interactions)
if setup.filter_features and setup.tech == "transcriptomics" and setup.data_type == "RSEM":

    # Transform to logCPM
    norm_factors = calculate_tmm_norm_factors(data.X)
    cpm_data = cpm(data.X, norm_factors = norm_factors, log = True)
    # Filter by signal/count
    min_counts = signal_by_percentile(cpm_data, percentile = 25)
    filtered_features = filter_by_signal(cpm_data, min_counts)
    # Subset
    filtered_data = data[:, filtered_features]

    # Variance filtering
    norm_factors = calculate_tmm_norm_factors(filtered_data.X)
    cpm_data = cpm(filtered_data.X, norm_factors = norm_factors, log = True)
    # Filter by variance
    min_variance = variance_by_percentile(cpm_data, percentile = 25)
    filtered_features = filter_by_variance(cpm_data, var_threshold = min_variance)

    # Subset features
    filtered_data = filtered_data[:, filtered_features]

elif setup.filter_features and setup.tech == "methylation" and setup.data_type == "Beta":

    # Transform to mvalues
    mvals = beta_to_mvalue(data.X)
    
    # Filter by variance
    min_variance = variance_by_percentile(mvals, percentile = 25)
    filtered_features = filter_by_variance(mvals, var_threshold = min_variance)

    # Subset features
    filtered_data = data[:, filtered_features]

elif setup.filter_features and setup.data_type == "proteomics":

    # No filtering
    filtered_data = data

# Save features (for DGE)
with open(setup.out("Features.yml"), "w") as f:
    yaml.dump(filtered_data.var_names.to_list(), f)

# Number of features
n_features = filtered_data.shape[1]
setup.add_settings({"n_features": n_features})

# Define training and inferred ancestry
setup.log("Define ancestries")
eur_data = filtered_data[filtered_data.obs[setup.ancestry_column] == setup.train_ancestry]
inf_data = filtered_data[filtered_data.obs[setup.ancestry_column] == setup.infer_ancestry]

# Assertion: Check min sample sizes per class in the eur data
data_validator.check_min_samples_per_class(
    data=eur_data.obs,
    column=setup.output_column,
    min_samples=setup.sample_cutoff,
    data_name="All train data"
    )


setup.log("Create subset")
# Classify covariates
covariate_list = setup.get("covariate", None)
if covariate_list:
    covariate_types = classify_covariates(eur_data.obs, covariates=covariate_list)
else:
    covariate_types = {"continuous": [], "discrete": []}

# Stratify based on output column and discrete covariates
stratification = [setup.output_column] + covariate_types["discrete"]
# Create dictionary with frequencies
strata = inf_data.obs.groupby(stratification, observed=False).size().to_dict()

# Stratified subset
train_data, test_data = stratified_subset(
    eur_data, 
    group_columns=stratification, 
    freq_dict=strata, 
    seed=seed
    )
train_idx = train_data.obs_names
test_idx = test_data.obs_names
inf_idx = inf_data.obs_names
# Assertion: Check for data leakage
data_validator.check_data_leakage(train_idx, test_idx)

# Assertion: Check min sample size per class
# Train data
data_validator.check_min_samples_per_class(
    data=train_data.obs,
    column=setup.output_column,
    min_samples=2,
    data_name="Train data"
    )
# Test data
data_validator.check_min_samples_per_class(
    data=test_data.obs,
    column=setup.output_column,
    min_samples=2,
    data_name="Test data"
    )
# Inf data
data_validator.check_min_samples_per_class(
    data=inf_data.obs,
    column=setup.output_column,
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

# Save observations
save_str = f"Obs_train.yml"
with open(setup.out(save_str), 'w') as f:
    yaml.dump(train_idx.to_list(), f)
# Add to settings
n_train = len(train_idx)
setup.add_settings({"n_train_obs": n_train})

save_str = f"Obs_test.yml"
with open(setup.out(save_str), 'w') as f:
    yaml.dump(test_idx.to_list(), f)
# Add to settings
n_test = len(test_idx)
setup.add_settings({"n_test_obs": n_test})

save_str = f"Obs_inf.yml"
with open(setup.out(save_str), 'w') as f:
    yaml.dump(inf_idx.to_list(), f)
# Add to settings
n_inf = len(inf_idx)
setup.add_settings({"n_inf_obs": n_inf})


# RUN R define R to run
R = ScriptRunner(r_path="/usr/bin/Rscript", py_path="/opt/conda/envs/ancestry/bin/python3")
settings_file = os.path.join(os.getcwd(), setup.output_directory, "Settings.yml")

# Statistical analysis (DGE)
setup.log("DGE")
script = "cross_ancestry_dge.R"
R.run_script(script_path=os.path.join(os.getcwd(), script), args=[settings_file])

# Machine learning
setup.log("Sklearn data")
# Encoding for the output column
comparison = setup.comparison
encoder = {j: i for i, j in enumerate(comparison)}
inv_encoder = {v: k for k, v in encoder.items()}

# Transform anndata to np.array
output_column = setup.output_column
train_X, train_y = np.array(train_data.X), encode_y(train_data.obs[output_column], encoder)
test_X, test_y = np.array(test_data.X), encode_y(test_data.obs[output_column], encoder)
inf_X, inf_y = np.array(inf_data.X), encode_y(inf_data.obs[output_column], encoder)

# Assertion: 
# Dimanesion of vectorsf
assert_str = "Dimension of feature matrix and true labes must align"
assert train_X.shape[0] == train_y.shape[0], assert_str
assert test_X.shape[0] == test_y.shape[0], assert_str
assert inf_X.shape[0] == inf_y.shape[0], assert_str
# Number of true labels
assert np.unique(train_y).shape[0] == 2
assert np.unique(test_y).shape[0] == 2
assert np.unique(inf_y).shape[0] == 2

# Normalization of features 
setup.log("Normalization")
# Select normalization method
tech = setup.tech
ml_normalization = setup.ml_normalization
normalization_method = ml_normalization_methods[tech][ml_normalization]
# Normalization
train_X = normalization_method(train_X)
test_X = normalization_method(test_X)
inf_X = normalization_method(inf_X)


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
        train_X, 
        train_y, 
        jobs_per_algorithm
        )
    for algo_name in setup.algorithms
)

# Evaluation
y_hat_list, interpretations_list, hyperparameters_list = [], [], []
for algo_name, best_model in zip(setup.algorithms, best_models):
    if best_model:
        y_hat, interpretations, hyperparameters = evaluate_algorithm_cross_ancestry(
            algo_name, 
            best_model,
            setup, 
            test_X, 
            test_y, 
            inf_X, 
            inf_y, 
            train_data.var_names,
            encoder = inv_encoder
        )
          
        # Append across algorithms
        y_hat_list.append(y_hat)
        interpretations_list.append(interpretations)
        hyperparameters_list.append(hyperparameters)
    else:
        print(f"Skipping evaluation for {algo_name} because training failed.")

# Combine across algorithms
# Probabilities
y_hat = pd.concat(y_hat_list)
# Add information
y_hat["Ancestry"] = setup.classification["infer_ancestry"]
y_hat["n_train_ancestry"] = len(train_idx)
y_hat["n_inf_ancestry"] = len(inf_idx)
y_hat["data_type"] = setup.data_type
# Save
y_hat.to_csv(setup.out(f"Probabilities.csv"), index=False)

# Interpretations
interpretaions = pd.concat(interpretations_list)
# Add information
interpretaions["Ancestry"] = setup.classification["infer_ancestry"]
interpretaions["n_train_ancestry"] = len(train_idx)
interpretaions["n_inf_ancestry"] = len(inf_idx)
interpretations["data_type"] = setup.data_type
# Save
interpretaions.to_csv(setup.out(f"Interpretations.csv"), index=False)

# Hyperparamters
hyperparameters = pd.concat(hyperparameters_list)
# Add information
hyperparameters["Ancestry"] = setup.classification["infer_ancestry"]
hyperparameters["n_train_ancestry"] = len(train_idx)
hyperparameters["n_inf_ancestry"] = len(inf_idx)
hyperparameters["data_type"] = setup.data_type
# Save
hyperparameters.to_csv(setup.out(f"Hyperparameters.csv"), index=False)

# Calculate metric
# Multiclass calculation
if setup.classification["multiclass"]:
    print('Score metrics for multiclass needs to be developed.')
# Binary calculation
else:
    metric_across_algo = []
    for algo_name in y_hat["Algorithm"].unique():
        # Filter by algorithm
        algo_probabilities = y_hat[y_hat["Algorithm"].str.contains(algo_name)]
        # Filter by status
        test_probabilities = algo_probabilities[algo_probabilities["Status"].str.contains("Test")]
        inf_probabilities = algo_probabilities[algo_probabilities["Status"].str.contains("Inference")]
        # Filter for target/class 1
        test_y, test_probabilities = test_probabilities.iloc[:,2].values, test_probabilities.iloc[:,1].values
        inf_y, inf_probabilities = inf_probabilities.iloc[:,2].values, inf_probabilities.iloc[:,1].values

        # ROC AUC
        test_auc_score = roc_auc_score(y_true=test_y, y_score=test_probabilities)
        inf_auc_score = roc_auc_score(y_true=inf_y, y_score=inf_probabilities)

        # Accuracy
        threshold = 0.5
        test_preds = (test_probabilities >= threshold).astype(int)
        inf_preds = (inf_probabilities >= threshold).astype(int)

        test_acc_score = accuracy_score(y_true=test_y, y_pred=test_preds)
        inf_acc_score = accuracy_score(y_true=inf_y, y_pred=inf_preds)

        # Prepare data
        train_data = {"Algorithm": algo_name, 
                      "ROC_AUC": test_auc_score,
                      "Accuracy": test_acc_score,
                      "Status": "Test", 
                      "Prediction": "Subset",
                      "Threshold": threshold
                    }
        inf_data = {"Algorithm": algo_name,
                    "ROC_AUC": inf_auc_score,
                    "Accuracy": inf_acc_score,
                    "Status": "Inference", 
                    "Prediction": "Ancestry",
                    "Threshold": threshold
                    }
        # Create metric dataframes
        test_metric = pd.DataFrame(train_data, index=[0])
        inf_metric = pd.DataFrame(inf_data, index=[0])
        # Combine
        metric_df = pd.concat([test_metric, inf_metric])
        # Add to across algorithms
        metric_across_algo.append(metric_df)
    
    # Join metric dataframes
    metric_across_algo = pd.concat(metric_across_algo)
    # Add additional information
    metric_across_algo["Seed"] = setup.seed
    metric_across_algo["Ancestry"] = setup["classification"]["infer_ancestry"].upper()
    metric_across_algo["n_train_ancestry"] = len(train_idx)
    metric_across_algo["n_inf_ancestry"] = len(inf_idx)
    metric_across_algo["data_type"] = setup.data_type
    # Save
    output_directory = setup["output_directory"]
    metric_across_algo.to_csv(os.path.join(output_directory, "Contrast_metric_ml.csv"), index=False)

# # Visualization
# if setup.visual_val:
#     setup.log('ML visualization')

#     # Run the script to visualize
#     R.run_script(script_path=os.path.join(os.getcwd(), "run_visualization.R"), 
#                  args=[settings_file])

# Save settings (overwrite with new settings)
with open(setup.out("Settings.yml"), 'w') as f:
    yaml.dump(setup.final_config, f)
setup.log("Settings done")


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