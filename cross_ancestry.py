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
    YAML_FILE = 'dev_settings.yml'
    print('Running interactive mode for development.')

# Here starts the script   

# Settings
setup = Setup(YAML_FILE)
with open(setup.out('Settings.yml'), 'w') as f:
    yaml.dump(setup.final_config, f)
setup.log('Settings done')


# Data setup
setup.log('Check data compatability')
data = anndata.read_h5ad(setup.data_path)


# Check if defined settings are present in the data
compatability = DataValidator(data=data, setup=setup)
compatability.validate()

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

# Assertion: Check for enough samples per class
# 1. Check for enough samples in the training ancestries (EUR)
# 2. Check for enough samples in inferred ancestry
counts = eur_data.obs[setup.classification['output_column']].value_counts()
          
assert_str = f'Prerequisite is a train sample size of {setup.sample_cutoff} per class.'
assert (counts > setup.sample_cutoff).all(), assert_str

# Check compared ancestry (infer_ancestry) if there are replicates
# Limma eBayes needs at least one replicate per class
counts = inf_data.obs[setup.classification['output_column']].value_counts()

assert_str = f'For DGE analysis of compared ancestry ({setup.classification['infer_ancestry'].upper()}) at least one replicate per class is needed.'
assert (counts >= 2).all(), assert_str

setup.log('Setting seed')
# Seed is used to generate different European subsets
seed = setup.seed
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

setup.log('Create subset')
# Proportions of output in inferred ancestry
freq = (inf_data.obs[setup.classification['output_column']]
        .value_counts()
        .to_dict()
        ) 
# Split the training data into train and test set (Based on frequency of compared ancestry)
train_data, test_data = stratified_subset(eur_data, freq, setup.classification['output_column'], seed)
train_idx, test_idx, inf_idx = train_data.obs_names, test_data.obs_names, inf_data.obs_names
# Assertion: Check for data leackage
assert_str = f'Data leakage! Observation also occur in the testset.'
assert not test_idx.isin(train_idx).all(), assert_str

# TODO - Check setup.covariate

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
# TODO - Use a non linear model e.g. Random Forest
for algo_name in setup.algorithms:

    setup.log(algo_name)
    # Load the sklearn class
    module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
    algo_class = getattr(importlib.import_module(module_name), algo_name)
    # Create instance of the algortihm
    algo_instance = algo_class(random_state=seed)


    # Cross-validation for hyperparameter search
    cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
    # Load the grid_search parameters
    search_space = setup.grid_search[algo_name]

    setup.log("Training")
    # Train with grid search and extract best model
    grid_search = GridSearchCV(estimator=algo_instance, param_grid=search_space,
                               cv=cv_hyp, scoring='f1_weighted', 
                               refit=True, n_jobs=setup.njobs
                               )
    grid_search.fit(train_X, train_y)
    best_m = grid_search.best_estimator_ 

    # Save hyperparameters 
    hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
    hyp.to_csv(setup.out(f"Hyperparameters_{algo_name}.csv"), index=False)

    # Save model
    if setup.save_model:
        with open(setup.out(f"{algo_name}_{setup.id}.pkl"), 'wb') as f:
            pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Training finished
    print(f"{algo_name} training done.")

    setup.log("Testing")
    # Validate prediction performance on EUR-subset
    test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
    test_y_hat["y"] = test_y
    # Save propabilities
    test_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_test.csv"), index=False)

    setup.log("Infering")
    # Validate prediction performance on ancestry
    inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
    inf_y_hat["y"] = inf_y
    # Save propabilities
    inf_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_inf.csv"), index=False)

    # Validation finished
    print(f"{algo_name} validation done")

    setup.log("Model interpretations")
    # Model feature importances
    feature_names = train_data.var_names
    feature_importance = extract_feature_importance(best_m, feature_names)
    # Save
    feature_importance.to_csv(setup.out(f"Feature_importance_{algo_name}.csv"), index=False)





# Algorithm (focusing on Logisitic Regression)
algo = LogisticRegression(random_state=seed)
# Cross-validation for hyperparameter search
cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
# Search space
search_space = setup.grid_search

setup.log('Start training')
grid_search = GridSearchCV(estimator=algo, param_grid=search_space,
                           cv=cv_hyp, scoring='f1_weighted', 
                           refit=True, n_jobs=setup.njobs
                           )
# Training the model
grid_search.fit(train_X, train_y)

# Best model
best_m = grid_search.best_estimator_ # is there a better version getter function

# Save hyperparameters 
hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
hyp.to_csv(setup.out('Hyperparameters.csv'), index=False)

# Save model
if setup.save_model:
    with open(setup.out(f'{setup.id}.pkl'), 'wb') as f:
        pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

# Training fisnished
print(f'Id: {setup.id} finished training.')

# Test the model on eur subset
setup.log('Testing')
test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
# Save propabilities
test_y_hat['y'] = test_y
test_y_hat.to_csv(setup.out('Probabilities_test.csv'), index=False)


# Test trained model on inferred ancestry 
setup.log('Inference')
# Trained model predicts inferred ancestry
inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
inf_y_hat['y'] = inf_y
# Save propabilities
inf_y_hat.to_csv(setup.out('Probabilities_inf.csv'), index=False)



# Model weights (interpretation)
feature_names = train_data.var_names # using train data because features got selected based on DE genes
feature_weights = best_m.coef_
feature_df = pd.DataFrame(data=feature_weights,
                          columns=feature_names)
feature_df.to_csv(setup.out('Weights.csv'), index=False)


# Calculate metric
# Multiclass calculation
if setup.classification['multiclass']:
    # TODO - Create multiclass scoring script
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

