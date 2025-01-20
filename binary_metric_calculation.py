# Standard libraries
import yaml
import pandas as pd
import numpy as np
import os
import sys
# Metric calculation
from sklearn.metrics import roc_auc_score

# Load settings
if len(sys.argv) > 1:
    file_path = sys.argv[1]
    with open(file_path, 'r') as f:
        setup = yaml.safe_load(f)
else:
    # Dev settings 
    file_path = 'dev_settings.yml'
    with open(file_path, 'r') as f:
        setup = yaml.safe_load(f)
    print('Running interactive mode for development.')

# Here starts the script
print('Start calculating metric.')

# Vscratch_dir is the place where the files are stored
vscratch_dir_in = setup["output_directory"]

# Number of samples
file_path = os.path.join(vscratch_dir_in, "Obs_train.yml")
with open(file_path, "r") as file:
    n_train_ancestry = len(yaml.safe_load(file))

file_path = os.path.join(vscratch_dir_in, "Obs_inf.yml")
with open(file_path, "r") as file:
    n_inf_ancestry = len(yaml.safe_load(file))

# Load the propabilities
metric_list = list()
for algo_name in setup["algorithms"]:

    # Names of dataframes
    test_df = f"Probabilities_{algo_name}_test.csv"
    inf_df = f"Probabilities_{algo_name}_inf.csv"

    # Load the propabilities
    test_probabilities = pd.read_csv(os.path.join(vscratch_dir_in, test_df))
    inf_probabilities = pd.read_csv(os.path.join(vscratch_dir_in, inf_df))

    # Filter for target/class 1
    test_y, test_probabilities = test_probabilities.iloc[:,-1].values, test_probabilities.iloc[:,1].values
    inf_y, inf_probabilities = inf_probabilities.iloc[:,-1].values, inf_probabilities.iloc[:,1].values

    # Calculate ROC AUC
    test_auc_score = roc_auc_score(y_true=test_y, y_score=test_probabilities)
    inf_auc_score = roc_auc_score(y_true=inf_y, y_score=inf_probabilities)
    
    # Prepare data
    train_data = {algo_name: test_auc_score, "Status": "Test", "Prediction": "Subset"}
    inf_data = {algo_name: inf_auc_score, "Status": "Inference", "Prediction": "Ancestry"}
    # Create metric dataframes
    test_metric = pd.DataFrame(train_data,index=[0])
    inf_metric = pd.DataFrame(inf_data,index=[0])
    # Combine
    metric_df = pd.concat([test_metric, inf_metric])

    # Add to list
    metric_list.append(metric_df)

# Join metric dataframes
merged_df = metric_list[0]
# Iterate over the remaining dataframes and merge them
for df in metric_list[1:]:
    merged_df = pd.merge(merged_df, df, on=["Status", "Prediction"], how="outer")

# Add additional information
merged_df = merged_df.assign(
    Seed=setup["seed"],
    Ancestry=setup["classification"]["infer_ancestry"].upper(),
    Metric="ROC_AUC",
    n_inf_ancestry=n_inf_ancestry,
    n_train_ancestry=n_train_ancestry
)

# Save
merged_df.to_csv(os.path.join(setup['output_directory'], 'Metric_ml.csv'), index=False)
