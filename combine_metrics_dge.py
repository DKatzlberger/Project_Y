# Standard libraries 
import os
import re
import pandas as pd
import numpy as np

# Command line parser
import argparse
import sys

# Metric
from sklearn.metrics import roc_auc_score

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
    YAML_FILE = 'job_settings.yml'
    print('Running interactive mode for development.')


with open(YAML_FILE, 'r') as f:
    data = yaml.safe_load(f)

# Modify output directory
base_dir = 'data/runs'
if data['tag'] == 'dev':
    directory = 'Dev'

    # Add path to the settings
    data['output_directory'] = os.path.join(base_dir, directory)
else:
    # Create a new direcory for each seed
    tag = data['tag']
    comparison = "_vs_".join(data['classification']['comparison'])
    comparison = f"{tag}_{comparison}"
    train_ancestry = data['classification']['train_ancestry'].upper()
    infer_ancestry = data['classification']['infer_ancestry'].upper()
    directory = f"{comparison}_{train_ancestry}_to_{infer_ancestry}_"

    # Add path to the settings
    data['output_directory'] = os.path.join(base_dir, directory)


train_data_frames = []
test_data_frames = []
inf_data_frames = []

# Loop through folders in the base directory
for folder_name in os.listdir(base_dir):
    # Check if the folder name follows the specific pattern and ends with an integer
    match = re.match(rf"{re.escape(directory)}(\d+)$", folder_name)
    if match:
        folder_path = os.path.join(base_dir, folder_name)
        if os.path.isdir(folder_path):
            # Extract the integer at the end of the folder name
            seed = int(match.group(1))

             # Load Contrast_train.csv if it exists
            train_file_path = os.path.join(folder_path, 'Contrast_train.csv')
            if os.path.exists(train_file_path):
                train_df = pd.read_csv(train_file_path)
                train_df['Status'] = 'Train'  # Add a column to identify the type
                train_df['Seed'] = seed  # Add the extracted folder number
                train_data_frames.append(train_df)

            test_file_path = os.path.join(folder_path, 'Contrast_test.csv')
            if os.path.exists(test_file_path):
                test_df = pd.read_csv(test_file_path)
                test_df['Status'] = 'Test'  # Add a column to identify the type
                test_df['Seed'] = seed  # Add the extracted folder number
                test_data_frames.append(test_df)

            inf_file_path = os.path.join(folder_path, 'Contrast_inf.csv')
            if os.path.exists(inf_file_path):
                inf_df = pd.read_csv(inf_file_path)
                inf_df['Status'] = 'Inf'  # Add a column to identify the type
                inf_df['Seed'] = seed  # Add the extracted folder number
                inf_data_frames.append(inf_df)


# Concatenate all data frames into one for each file type
all_train_data = pd.concat(train_data_frames, ignore_index=True) 
all_test_data = pd.concat(test_data_frames, ignore_index=True) 
all_inf_data = pd.concat(inf_data_frames, ignore_index=True) 

# Merge dataframes for correlation
train_data = all_train_data[['coef', 'Gene', 'logFC', 'Seed']]
test_data = all_test_data[['coef', 'Gene', 'logFC', 'Seed']]
inf_data = all_inf_data[['coef', 'Gene', 'logFC', 'Seed']]

merged_train_test = train_data.merge(test_data, 
                                     on=['Gene', 'coef', 'Seed'], 
                                     suffixes=('_train', '_test')
                                     )

merged_data = (merged_train_test
               .merge(inf_data, on=['Gene', 'coef', 'Seed'])
               .rename(columns={'logFC': 'logFC_inf'})
               .drop(columns=['Seed'])
               )

# Correlate logFCs
pearson_corr = merged_data.corr(numeric_only=True, method='pearson')
spearman_corr = merged_data.corr(numeric_only=True, method='spearman')

# Rename index
pearson_corr.index = pearson_corr.index.str.replace("logFC_", "")
pearson_corr.columns = pearson_corr.columns.str.replace("logFC_", "")

# Melt the matrix
upper_triangle = pd.DataFrame(pearson_corr.where(np.triu(np.ones(pearson_corr.shape), k=1).astype(bool)).iloc[0].dropna())
df = upper_triangle.reset_index().rename(columns={'index': 'Status',
                                                  'train': 'Pearson_correlation'})

df = df.assign(
    **{'Ancestry': infer_ancestry,
       'Comparison': comparison})
df['Prediction'] = np.where(df['Status'] == 'test', 'Subset', 'Ancestry')

# Save 
base_output_dir = 'data/combined_runs'
extension = 'combined_dge_metric'
output_csv_path = os.path.join(base_output_dir, directory + extension + '.csv')
df.to_csv(output_csv_path, index=False)