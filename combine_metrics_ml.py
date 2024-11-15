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


# Machine learning 
# Initialize lists to collect data frames
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
            
            # Load Probabilities_test.csv if it exists
            test_file_path = os.path.join(folder_path, "Probabilities_test.csv")
            if os.path.exists(test_file_path):
                test_df = pd.read_csv(test_file_path)
                test_df['Seed'] = seed  # Add metadata column
                test_df['Status'] = 'Test'     # Label the file type
                test_data_frames.append(test_df)  # Collect the DataFrame
                

            # Load Probabilities_inf.csv if it exists
            inf_file_path = os.path.join(folder_path, "Probabilities_inf.csv")
            if os.path.exists(inf_file_path):
                inf_df = pd.read_csv(inf_file_path)
                inf_df['Seed'] = seed  # Add metadata column
                inf_df['Status'] = 'Inf'      # Label the file type
                inf_data_frames.append(inf_df)   # Collect the DataFrame
                

# Concatenate all data frames into one for each file type
all_test_data = pd.concat(test_data_frames, ignore_index=True) 
all_inf_data = pd.concat(inf_data_frames, ignore_index=True) 

# Calculate ROC AUC
auc_scores = {}
test_auc = roc_auc_score(all_test_data.iloc[:,2], all_test_data.iloc[:,1])
auc_scores['Test'] = test_auc

inf_auc = roc_auc_score(all_inf_data.iloc[:,2], all_inf_data.iloc[:,1])
auc_scores['Inf'] = inf_auc

# Create the DataFrame as per your structure
df = (pd.DataFrame.from_dict(auc_scores, orient='index', columns=['ROC_AUC'])
      .assign(
          **{'Ancestry': infer_ancestry,
             'Comparison': comparison}
          )
      .reset_index(names='Status')
)
df['Prediction'] = np.where(df['Status'] == 'Test', 'Subset', 'Ancestry')
df['Analysis'] = 'Machine_learning'

# Save the file
base_output_dir = 'data/combined_runs'
extension = 'combined_ml_metric'
output_csv_path = os.path.join(base_output_dir, directory + extension + '.csv')
df.to_csv(output_csv_path, index=False)
