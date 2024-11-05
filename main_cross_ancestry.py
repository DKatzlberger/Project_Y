# Standard libraries
import yaml
import uuid
import pandas as pd
import numpy as np
import anndata
import pickle

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

# Here starts the script   

# Settings
setup = Setup('dev_settings.yml')
with open(setup.out('Settings.yml'), 'w') as f:
    yaml.dump(setup.final_config, f)
setup.log('Settings done')

# Data setup
data = anndata.read_h5ad(setup.data_path)

# Make feature names unique
# data.var_names_make_unique()
assert_str = 'Features not unique.'
assert data.var.index.shape[0] == len(set(data.var.index)), assert_str

# Define training and inferred ancestry
setup.log('Define ancestries')
eur_data = data[data.obs.genetic_ancestry == setup.classification['train_ancestry']]
inf_data = data[data.obs.genetic_ancestry == setup.classification['infer_ancestry']]
# Console output
# print(f'Algorithms are trained on {setup.classification['train_ancestry'].upper()}' 
#       f' and inferred on {setup.classification['infer_ancestry'].upper()}.')

# Define classification task 
setup.log('Define classification')
eur_data = (eur_data[eur_data.obs[setup.classification['output_column']]
                     .isin(setup.classification['comparison'])]
                     )
inf_data = (inf_data[inf_data.obs[setup.classification['output_column']]
                     .isin(setup.classification['comparison'])]
                     )

# Assertion: Check for enough samples per class
counts = (eur_data.obs[setup.classification['output_column']]
          .value_counts()
          ) 
assert_str = f'Prerequisite is a sample size of {setup.sample_cutoff} for each class.'
assert (counts > setup.sample_cutoff).all(), assert_str

# Normalization of features 
# setup.log('Normalization')
# eur_data.X = normalize(eval(setup.normalization), eur_data.X)
# inf_data.X = normalize(eval(setup.normalization), inf_data.X) 

setup.log('Setting seed')
# Seed is used to generate different European subsets
seed = setup.seed
np.random.seed(seed)
# os.environ['PYTHONHASHSEED'] = str(seed)
# check this sklearn seed 

setup.log('Data setup')
# Proportions of output in inferred ancestry
freq = (inf_data.obs[setup.classification['output_column']]
        .value_counts()
        .to_dict()
        ) 
# Split the training data into train and test set (Based on frequency of compared ancestry)
train_data, test_data = stratified_subset(eur_data, freq, setup.classification['output_column'], seed)
train_idx, test_idx, inf_idx = train_data.obs_names, test_data.obs_names, inf_data.obs_names
# Assertion: Check for data leackage
assert_str = f'Observation also occur in the testset.'
assert not test_idx.isin(train_idx).all(), assert_str

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

# Statistical analysis (DGE)
setup.log('DGE')
r_input = os.path.join(setup.output_directory, 'Settings.yml')
r_execute = '/usr/bin/Rscript'
path2script = os.path.join(os.getcwd(), 'dge_workflow.R')
# Calling R script to do analysis
subprocess.run([r_execute, path2script, r_input])

# Subset features that are used in DGE
p = os.path.join(setup.output_directory, 'Features.yml')
dge_features = yaml.safe_load(Path(p).read_text())

train_data = train_data[:, dge_features]
test_data = test_data[:, dge_features]
inf_data = inf_data[:, dge_features]

# Machine learning
print('Starting machine learning.')
setup.log('Sklearn data')
# Transform anndata to np.array
train_X = np.array(train_data.X)
train_y, _ = create_y(train_data, setup.classification['output_column'])

test_X = np.array(test_data.X)
test_y, _ = create_y(test_data, setup.classification['output_column']) 

# Assertion: Check arrays
# Labels
assert_str = f'There have to be at least 2 labels, {np.unique(train_y).shape[0]} are given.'
assert np.unique(train_y).shape[0] >= 2, assert_str
# Dimension
assert_str = f'Dimensions of label vector ({train_y.shape[0]}) and feature matrix ({train_X.shape[0]}) do not align.'
assert train_y.shape[0] == train_X.shape[0], assert_str

setup.log('Model setup')
# Algorithm (focusing on Logisitic Regression)
algo = LogisticRegression(random_state=seed)
# Cross-validation for hyperparameter search
cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=seed)
# Search space
search_space = setup.grid_search

setup.log('Start training')
grid_search = GridSearchCV(estimator=algo, param_grid=search_space,
                           cv=cv_hyp, scoring='f1_weighted', 
                           refit=True
                           )
# Training the model
grid_search.fit(train_X, train_y)

# Best model
best_m = grid_search.best_estimator_ # is there a better version getter function

# Test the model on the eur subset
test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
# Save propabilities
test_y_hat['y'] = test_y
test_y_hat['ancestry'] = setup.classification['train_ancestry'].upper()
test_y_hat.to_csv(setup.out('Probabilities_test.csv'), index=False)

# Save hyperparameters 
hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
hyp.to_csv(setup.out('Hyperparameters.csv'), index=False)

# Save model
if setup.save_model:
    with open(setup.out(f'{setup.id}.pkl'), 'wb') as f:
        pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL) # check if there is a better function

# Training fisnished
print(f'{setup.id} finished')

# Test trained model on different ancestry to check if learned patterns generalize
setup.log('Inference')
inf_X = inf_data.X
inf_y, _ = create_y(inf_data, setup.classification['output_column'])
# Trained model predicts inferred ancestry
inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
inf_y_hat['y'] = inf_y
inf_y_hat['ancestry'] = setup.classification['infer_ancestry'].upper()
# Save propabilities
inf_y_hat.to_csv(setup.out('Probabilities_inf.csv'), index=False)

# Plot for visualization
# TODO - Outsource to a R file
# Calculate ROC AUC as metric
test_auc = roc_auc_score(test_y, best_m.predict_proba(test_X)[:, 1])
inf_auc = roc_auc_score(inf_y, best_m.predict_proba(inf_X)[:, 1])

# Dataframe for plotting
df = pd.concat(
    [pd.DataFrame({'ROC AUC': [test_auc], 
                   'Prediction': ['European subset']
                   }, 
                   index=[0]
                   ), 
     pd.DataFrame({'ROC AUC': [inf_auc],
                   'Prediction': [f'Ancestry\n({setup.classification['infer_ancestry'].upper()})']
                   }, 
                   index=[1]
                   )]
                   )
# Plot
g = sns.barplot(data=df,
                x='Prediction',
                y='ROC AUC')
g.figure.savefig(f'{setup.out('Metric.pdf')}')

# TODO - Model interpretations