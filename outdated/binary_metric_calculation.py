# Standard libraries
import yaml
import pandas as pd
import numpy as np
import os
import sys
# Metric calculation
from sklearn.metrics import roc_auc_score, accuracy_score

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

# Load the propabilities
file = "Probabilities.csv"
probabilities = pd.read_csv(os.path.join(vscratch_dir_in, file))

# Number of observation
n_train_ancestry = probabilities["n_train_ancestry"]
n_inf_ancestry = probabilities["n_inf_ancestry"]

metric_list = list()
for algo_name in probabilities["Algorithm"].unique():

    # Filter by algorithm
    algo_probabilities = probabilities[probabilities["Algorithm"].str.contains(algo_name)]

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

    # Add to list
    metric_list.append(metric_df)

# Join metric dataframes
merged_df = pd.concat(metric_list)

# Add additional information
merged_df = merged_df.assign(
    Seed=setup["seed"],
    Ancestry=setup["classification"]["infer_ancestry"].upper(),
    n_inf_ancestry=n_inf_ancestry,
    n_train_ancestry=n_train_ancestry
)

# Save
merged_df.to_csv(os.path.join(setup['output_directory'], 'Contrast_metric_ml.csv'), index=False)
