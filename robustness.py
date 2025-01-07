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

# Subprocess to run R
import subprocess

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


# Load the settings
# Don't need setup class?
setup = Setup(YAML_FILE)
with open(setup.out('Settings.yml'), 'w') as f:
    yaml.dump(setup.final_config, f)
setup.log('Settings done')

# Load the data
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

# Validate classification task
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

    # Use class frequency to mimick inferred ancestry
    freq = inf_downsampled.obs[setup.classification['output_column']].value_counts().to_dict()
    # Create 'train_data' and 'test_data' (EUR-subset)
    train_data, test_data = stratified_subset(eur_downsampled, freq, setup.classification['output_column'], seed)
    # Get observations 
    train_idx, test_idx, inf_idx = train_data.obs_names, test_data.obs_names, inf_downsampled.obs_names

    # Assertion: 
    # Check for data leackage
    assert_str = f'Data leakage! Observation also occur in the testset.'
    assert not test_idx.isin(train_idx).all(), assert_str

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
robustness_results = list()
for prop in setup.proportion:
    # Filter data based on observations and features per proportion
    prop = str(prop).replace(".", "_")

    # Create file names for observations
    obs_file_train = os.path.join(setup.output_directory,
                                  f"Obs_proportion_{prop}_train.yml")
    obs_file_test = os.path.join(setup.output_directory,
                                 f"Obs_proportion_{prop}_test.yml")
    obs_file_inf = os.path.join(setup.output_directory,
                                f"Obs_proportion_{prop}_inf.yml")
    # Load observations
    obs_train = yaml.safe_load(Path(obs_file_train).read_text())
    obs_test = yaml.safe_load(Path(obs_file_test).read_text())
    obs_inf = yaml.safe_load(Path(obs_file_inf).read_text())

    # Feature file name
    feature_file = os.path.join(setup.output_directory,
                                f"Features_{prop}.yml")
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

    # Normalization of features 
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
    # Algorithm (focusing on Logisitic Regression)
    algo = LogisticRegression(random_state=seed)
    # Cross-validation for hyperparameter search
    cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
    # Search space
    search_space = setup.grid_search

    # Training
    grid_search = GridSearchCV(estimator=algo, param_grid=search_space,
                           cv=cv_hyp, scoring='f1_weighted', 
                           refit=True, n_jobs=setup.njobs
                           )
    # Training the model
    grid_search.fit(train_X, train_y)
    # Best model
    best_m = grid_search.best_estimator_ 

    # Save hyperparameters 
    hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
    hyp.to_csv(setup.out(f'Hyperparameters_{prop}.csv'), index=False)

    # Save model
    if setup.save_model:
        with open(setup.out(f'{setup.id}.pkl'), 'wb') as f:
            pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Validation of the model
    # Test the model on eur subset
    test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
    # Save propabilities
    test_y_hat['y'] = test_y
    test_y_hat.to_csv(setup.out(f'Probabilities_{prop}_test.csv'), index=False)

    # Test trained model on inferred ancestry 
    # Trained model predicts inferred ancestry
    inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
    inf_y_hat['y'] = inf_y
    # Save propabilities
    inf_y_hat.to_csv(setup.out(f'Probabilities_{prop}_inf.csv'), index=False)

    # Model weights (interpretation)
    feature_names = train_data.var_names # using train data because features got selected based on DE genes
    feature_weights = best_m.coef_
    feature_df = pd.DataFrame(data=feature_weights,
                            columns=feature_names)
    feature_df.to_csv(setup.out(f'Weights_{prop}.csv'), index=False)

    # Calculate metric
    # Test-set
    props, labels = np.array(test_y_hat.iloc[:, 1]), np.array(test_y_hat['y'])
    test_auc = roc_auc_score(y_true=labels, y_score=props)

    # Inf-set
    props, labels = np.array(inf_y_hat.iloc[:, 1]), np.array(inf_y_hat['y'])
    inf_auc = roc_auc_score(y_true=labels, y_score=props)

    # Metric dataframe
    test_metric = pd.DataFrame({"ROC_AUC": test_auc,
                                "Status": "Test",
                                "Prediction": "Subset"
                                },
                                index=[0])
    
    inf_metric = pd.DataFrame({"ROC_AUC": inf_auc,
                               "Status": "Inference",
                               "Prediction": "Ancestry"
                               },
                               index=[0])
    # Combine
    metric_df = (pd.concat([test_metric, inf_metric])
                 .assign(
                       **{"n_train_ancestry": len(obs_inf),
                          "n_test_ancetry": len(obs_train),
                          "proportion": prop.replace("_", "."),
                          "Ancestry": setup.classification["infer_ancestry"].upper()
                          }
                 ))

    robustness_results.append(metric_df)

# Combine all dataframes and save 
final_metric_df = pd.concat(robustness_results)
final_metric_df.to_csv(os.path.join(setup['output_directory'], 'Metric_ml.csv'), index=False)



