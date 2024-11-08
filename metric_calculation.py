# Standard libraries
import yaml
import pandas as pd
import numpy as np
import os
import sys
# Metric calculation
import sklearn.metrics as mc

# Functions 
def binary_metric_score(y_true, y_prob, thresholds=None, metric='accuracy'):
     
    metric_functions = {
        'accuracy': mc.accuracy_score,
        'f1': mc.f1_score,
        'roc_auc': mc.roc_auc_score
    }

    # Check if metric is valid
    if metric not in metric_functions:
        raise ValueError("Invalid metric specified. Choose 'accuracy', 'f1' or ''roc_auc.")
    
    # Special handling for ROC AUC
    if metric == 'roc_auc':
        return mc.roc_auc_score(y_true, y_prob)
    
    # Check if thresholds are given
    if thresholds is None:
        raise ValueError('Thresholds must be provided for accuracy and F1 score calculations.')

    # Select metric
    func = metric_functions[metric]

    # Dict for  scores
    scores = {}
    for thr in thresholds:
        # Convert to class labels
        y_hat = (y_prob >= thr).astype(int)

        # Calculate metric
        metric = func(y_true, y_hat)
        scores[thr] = metric
    
    return scores


# TODO - This is only availabe for binary classification

# Load settings
file_path = sys.argv[1]
with open(file_path, 'r') as f:
    setup = yaml.safe_load(f)

# Test probabilities (European subset)
test_p = pd.read_csv(os.path.join(setup['output_directory'], 'Probabilities_test.csv'))
# Take probabilities for class 1
test_y, test_p = np.array(test_p.iloc[:,-1]), np.array(test_p.iloc[:,1])

# Inf probabilities
inf_p = pd.read_csv(os.path.join(setup['output_directory'], 'Probabilities_inf.csv'))
# Take probabilities for class 1
inf_y, inf_p = np.array(inf_p.iloc[:,-1]), np.array(inf_p.iloc[:,1])

# Thresholds for dependent metric
thr = [0.2, 0.3, 0.4, 0.5]

# Metric calculation
test_auc = binary_metric_score(test_y, test_p, metric='roc_auc')
test_acc = binary_metric_score(test_y, test_p, thresholds=thr, metric='accuracy')
test_f1 = binary_metric_score(test_y, test_p, thresholds=thr, metric='f1')

inf_auc = binary_metric_score(inf_y, inf_p, metric='roc_auc')
inf_acc = binary_metric_score(inf_y, inf_p, thresholds=thr, metric='accuracy')
inf_f1 = binary_metric_score(inf_y, inf_p, thresholds=thr, metric='f1')

# Make dataframe
test_metric = (pd.concat(
    [pd.DataFrame.from_dict(metric, orient='index', columns=[name])
     for metric, name in zip([test_acc, test_f1], ['Accuracy', 'F1'])],
    axis=1
)
.reset_index(names='Threshold')
.assign(
    **{
    'ROC AUC': test_auc,
    'Prediction': 'Eur Subset',
    'Status': 'Test',
    'Ancestry': setup['classification']['infer_ancestry'].upper()
}))

inf_metric = (pd.concat(
    [pd.DataFrame.from_dict(metric, orient='index', columns=[name])
     for metric, name in zip([inf_acc, inf_f1], ['Accuracy', 'F1'])],
    axis=1
)
.reset_index(names='Threshold')
.assign(
    **{
    'ROC AUC': inf_auc,
    'Prediction': 'Ancestry',
    'Status': 'Inference',
    'Ancestry': setup['classification']['infer_ancestry'].upper()
}))

# This variable is then used in R
r_metric = pd.concat([test_metric, inf_metric]).reset_index(drop=True)

