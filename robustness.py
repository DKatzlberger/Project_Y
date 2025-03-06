# Currently the only difference to 'main_cross_ancestry'
# is the downsampling of the data to a certain proportion
# e.g. 90%

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
# Log
setup.log("Settings done")

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
setup.log("Define ancestries")
eur_data = data[data.obs[setup.classification['ancestry_column']] == setup.classification['train_ancestry']]
inf_data = data[data.obs[setup.classification['ancestry_column']] == setup.classification['infer_ancestry']]

# Define classification task 
setup.log("Define classification")
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

setup.log("Setting seed")
# Seed is used to generate different European subsets
# Subsets are sampled randomly
seed = setup.seed
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# To check robustness only use 90% or 80% of data 
# Train data is set based on the inferred ancestry
# 1. E.g. only 90% of data is used
#   a) Downsample to 0.9 of data
#   b) Based on inferred ancestry the EUR-subset is created 

# Set proportion to 0.9
# To downsample original data
proportions = setup.proportion

# Downsample 'eur_data' by proportion
setup.log('Downsample data')
eur_downsampled_list = sample_by_size(adata=eur_data,
                                      props=proportions,
                                      output_column=setup.classification['output_column'],
                                      seed=seed)

# Downsample 'inf_data' by proportion
inf_downsampled_list = sample_by_size(adata=inf_data,
                                      props=proportions,
                                      output_column=setup.classification['output_column'],
                                      seed=seed)

# Iterate over all downsampled datasets
# 1. Create EUR-subset based on 'inf_downsampeld'
# 2. Get observations and save them 
for downsampled_data in zip(eur_downsampled_list, inf_downsampled_list):
    eur_downsampled = eur_downsampled_list[downsampled_data[0]]
    inf_downsampled = inf_downsampled_list[downsampled_data[1]]

    # TODO - Covariant stuff
    covariate_list = setup.classification.get("covariate", None)
    if covariate_list:
        covariate_types = classify_covariates(eur_data.obs, covariates=covariate_list)
    else:
        covariate_types = {'continuous': [], 'discrete': []}
    
    # Stratify based on output column and discrete covariates
    to_stratify_columns = [setup.classification["output_column"]] + covariate_types["discrete"]
    # Create dictionary with frequencies
    frequencies = inf_downsampled.obs.groupby(to_stratify_columns, observed=False).size().to_dict()

    # Create 'train_data' and 'test_data' (EUR-subset)
    train_data, test_data = stratified_subset(eur_downsampled, group_columns=to_stratify_columns, freq_dict=frequencies, seed=seed)
    # Get observations 
    train_idx, test_idx, inf_idx = train_data.obs_names, test_data.obs_names, inf_downsampled.obs_names

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
        data=inf_downsampled.obs,
        column=setup.classification['output_column'],
        min_samples=2,
        data_name="Inf data"
        )

    # Check for correct downsampling
    assert_str = f"Train and test set do not add up to complete data."
    assert len(train_idx) + len(test_idx) == len(eur_downsampled), assert_str 

    # Save the observations
    # Remaining Europeans
    save_str = f'Obs_{downsampled_data[0]}_train.yml'
    with open(setup.out(save_str), 'w') as f:
        yaml.dump(train_idx.to_list(), f)

    # European subset
    save_str = f'Obs_{downsampled_data[0]}_test.yml'
    with open(setup.out(save_str), 'w') as f:
        yaml.dump(test_idx.to_list(), f)
    
    # Inferred Ancestry
    save_str = f'Obs_{downsampled_data[1]}_inf.yml'
    with open(setup.out(save_str), 'w') as f:
        yaml.dump(inf_idx.to_list(), f)


# Statistical analysis (DGE)
setup.log('DGE')
# ScriptRunner defines bin for R and python
R = ScriptRunner(r_path='/usr/bin/Rscript', py_path='/opt/conda/envs/ancestry/bin/python3')
settings_file = os.path.join(os.getcwd(), setup.output_directory, 'Settings.yml')

# Runs R script
script = "robustness_dge.R"
R.run_script(script_path=os.path.join(os.getcwd(), script),
             args=[settings_file])

# Iterate over all proportions and filter data
setup.log(f"Machine learning")
robustness_probabilities = list()
robustness_metric = list()
for prop in setup.proportion:
    # Filter data based on observations and features per proportion
    prop_ = str(prop).replace(".", "_")

    # Create file names for observations
    obs_file_train = os.path.join(setup.output_directory,
                                  f"Obs_proportion_{prop_}_train.yml")
    obs_file_test = os.path.join(setup.output_directory,
                                 f"Obs_proportion_{prop_}_test.yml")
    obs_file_inf = os.path.join(setup.output_directory,
                                f"Obs_proportion_{prop_}_inf.yml")
    # Load observations
    obs_train = yaml.safe_load(Path(obs_file_train).read_text())
    obs_test = yaml.safe_load(Path(obs_file_test).read_text())
    obs_inf = yaml.safe_load(Path(obs_file_inf).read_text())

    # Feature file name
    feature_file = os.path.join(setup.output_directory,f"Features_{prop_}.yml")
    # Load features
    dge_features = yaml.safe_load(Path(feature_file).read_text())

    # Filter data
    train_data = data[obs_train, dge_features]
    test_data = data[obs_test, dge_features]
    inf_data = data[obs_inf, dge_features]

    # Machine learning
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

    # Normalization
    # Select normalization method
    data_type = setup.data_type
    ml_normalization = setup.ml_normalization
    normalization_method = ml_normalization_methods[data_type][ml_normalization]
    # Normalization
    train_X = normalization_method(train_X)
    test_X = normalization_method(test_X)
    inf_X = normalization_method(inf_X)

    # Assertion: Check arrays
    # Labels
    assert_str = f'There have to be at least 2 classes, {np.unique(train_y).shape[0]} are given.'
    assert np.unique(train_y).shape[0] >= 2, assert_str
    # Dimension
    assert_str = f'Dimensions of label vector ({train_y.shape[0]}) and feature matrix ({train_X.shape[0]}) do not align.'
    assert train_y.shape[0] == train_X.shape[0], assert_str

    print(f"Start machine learning proportion: {prop}")
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
            jobs_per_algorithm, 
            prop=prop_
            )
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
                encoder = inv_encoder,
                prop = prop_
                )
            # Append probabilities
            y_hat_list.append(y_hat)
        else:
            print(f"Skipping evaluation for {algo_name} because training failed.")

    # Combine probabilities across algorithms
    y_hat = pd.concat(y_hat_list)
    # Add information
    y_hat["Proportion"] = prop
    y_hat["n_train_ancestry"] = len(train_idx)
    y_hat["n_inf_ancestry"] = len(inf_idx)
    # Append across proportions
    robustness_probabilities.append(y_hat)

    # Calculate metric
    metric_list = list()
    for algo_name in y_hat["Algorithm"].unique():
        # Filter by algorithm
        algo_probabilities = y_hat[y_hat["Algorithm"].str.contains(algo_name)]
        # Filter by status
        test_probabilities = algo_probabilities[algo_probabilities["Status"].str.contains("Test")]
        inf_probabilities = algo_probabilities[algo_probabilities["Status"].str.contains("Inference")]
        # Filter for target/class 1
        test_y, test_probabilities = test_probabilities.iloc[:,2].values, test_probabilities.iloc[:,1].values
        inf_y, inf_probabilities = inf_probabilities.iloc[:,2].values, inf_probabilities.iloc[:,1].values

        # Calculate ROC AUC
        test_auc_score = roc_auc_score(y_true=test_y, y_score=test_probabilities)
        inf_auc_score = roc_auc_score(y_true=inf_y, y_score=inf_probabilities)

        # Calculate Accuracy
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
        # Add across algorithms
        metric_list.append(metric_df)

    # Join metric dataframes
    merged_df = pd.concat(metric_list)
    # Add additional information
    merged_df["Proportion"] = prop
    merged_df["Seed"] = setup.seed
    merged_df["Ancestry"] = setup["classification"]["infer_ancestry"].upper()
    merged_df["n_train_ancestry"] = len(train_idx)
    merged_df["n_inf_ancestry"] = len(inf_idx)

    # Add to across proportions
    robustness_metric.append(merged_df)

# Combine across proportions
robustness_metric_ml = pd.concat(robustness_metric)
robustness_probabilities = pd.concat(robustness_probabilities)
# Save
output_directory = setup["output_directory"]
robustness_metric_ml.to_csv(os.path.join(output_directory, "Contrast_metric_ml.csv"), index=False)
robustness_probabilities.to_csv(os.path.join(output_directory, "Probabilities.csv"), index=False)




#     # Loop over algorithms
#     algo_list = list()
#     for algo_name in setup.algorithms:

#         setup.log(algo_name)
#         # Load the sklearn class
#         module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
#         algo_class = getattr(importlib.import_module(module_name), algo_name)
#         # Create instance of the algortihm
#         algo_instance = algo_class(random_state=seed)

#         # Cross-validation for hyperparameter search
#         cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
#         # Load the grid_search parameters
#         search_space = setup.grid_search[algo_name]

#         setup.log("Training")
#         # Train with grid search and extract best model
#         grid_search = GridSearchCV(estimator=algo_instance, param_grid=search_space,
#                                    cv=cv_hyp, scoring='f1_weighted', 
#                                    refit=True, n_jobs=setup.njobs
#                                     )
#         grid_search.fit(train_X, train_y)
#         best_m = grid_search.best_estimator_ 

#         # Save hyperparameters 
#         hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
#         hyp.to_csv(setup.out(f"Hyperparameters_{prop_}_{algo_name}.csv"), index=False)

#         # Save model
#         if setup.save_model:
#             with open(setup.out(f"{algo_name}_{prop_}_{setup.id}.pkl"), 'wb') as f:
#                 pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

#         # Training finished
#         print(f"{algo_name} training done.")

#         setup.log("Testing")
#         # Validate prediction performance on EUR-subset
#         test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
#         test_y_hat["y"] = test_y
#         # Save propabilities
#         test_y_hat.to_csv(setup.out(f"Probabilities_{prop_}_{algo_name}_test.csv"), index=False)

#         setup.log("Infering")
#         # Validate prediction performance on ancestry
#         inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
#         inf_y_hat["y"] = inf_y
#         # Save propabilities
#         inf_y_hat.to_csv(setup.out(f"Probabilities_{prop_}_{algo_name}_inf.csv"), index=False)

#         # Validation finished
#         print(f"{algo_name} validation done")

#         setup.log("Model interpretations")
#         # Extract feature importances
#         feature_names = train_data.var_names
#         feature_importance = extract_feature_importance(best_m, feature_names)
#         # Save feature importances
#         feature_importance.to_csv(setup.out(f"Feature_importance_{prop_}_{algo_name}.csv"), index=False)

#         # Calculate metric
#         test_y, test_probabilities = test_y_hat.iloc[:, -1].values, test_y_hat.iloc[:,1].values
#         inf_y, inf_probabilities = inf_y_hat.iloc[:,-1].values, inf_y_hat.iloc[:,1].values

#         # Calculate ROC AUC
#         test_auc_score = roc_auc_score(y_true=test_y, y_score=test_probabilities)
#         inf_auc_score = roc_auc_score(y_true=inf_y, y_score=inf_probabilities)

#         # Prepare data
#         train_data_metric = {algo_name: test_auc_score, "Status": "Test","Prediction": "Subset"}
#         inf_data_metric = {algo_name: inf_auc_score, "Status": "Inference","Prediction": "Ancestry"}
#         # Create metric dataframes
#         test_metric = pd.DataFrame(train_data_metric,index=[0])
#         inf_metric = pd.DataFrame(inf_data_metric,index=[0])
#         # Combine
#         metric_df = pd.concat([test_metric, inf_metric])

#         # Add to list
#         algo_list.append(metric_df)

#     # Join algorithm metric dataframes
#     merged_df = algo_list[0]
#     # Iterate over the remaining dataframes and merge them
#     for df in algo_list[1:]:
#         merged_df = pd.merge(merged_df, df, on=["Status", "Prediction"], how="outer")

#     # Add additional information 
#     merged_df = merged_df.assign(
#         Seed=seed,
#         n_inf_ancestry=len(obs_inf),
#         n_train_ancestry=len(obs_train),
#         Proportion=prop,
#         Ancestry=setup.classification["infer_ancestry"].upper()
#     )

#     # Add to list 
#     robustness_list.append(merged_df)

# # Combine all proportions
# robustness_results = pd.concat(robustness_list)
# # Save
# robustness_results.to_csv(setup.out("Metric_ml.csv"), index=False)