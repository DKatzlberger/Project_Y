# Standard libraries
import numpy as np
import pandas as pd
import yaml
import uuid
import os
import time
from pathlib import Path
from datetime import datetime
import psutil
import subprocess

# Object that contains settings and methods
class Setup(dict):
    def __init__(self, config_file):
        
        # TODO - Maybe combination of default and custom settings outsource
        # Open the default config file
        try:
            default_config = yaml.safe_load(Path('default_settings.yml').read_text())
        except yaml.YAMLError as exc:
            print(exc)

        # Open custom config file
        try:
            custom_config = yaml.safe_load(Path(config_file).read_text())
        except yaml.YAMLError as exc:
            print(exc)
        
        # Combine custom and default settings 
        final_config = {**default_config, **custom_config} 

        # add an id and date
        final_config['id'] = uuid.uuid4().hex.upper()[:10]
        final_config['date'] = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())

        # Print statement to inform user about their analysis
        print(f'New analysis with id: {final_config['id']}; created: {final_config['date']}')

        # Check required settings
        required_settings = ['classification', 'data_path', 'output_directory', 'seed']
        for i in required_settings:
            assert i in final_config, f'No {i} defined but required'     
        
        # Check required settings in classification
        required_settings = ['comparison', 'output_column', 'train_ancestry', 'infer_ancestry']
        for i in required_settings:
             assert i in final_config['classification'], f'No {i} in classification but required'


        # TODO - Check for multiclass
        # if len(final_config['classification']['comparison']) > 2:
        #     final_config['classification'].set      

        # Init
        self.final_config = final_config
        super().__init__(final_config)

        # Create output directory
        Path(final_config['output_directory']).mkdir(parents=True, exist_ok=self.overwrite)
        save_location = os.path.join(os.getcwd(), final_config['output_directory'])
        # Print tatement to inform about save location
        print(f'Output will be saved to {save_location}')

        self.check_settings()
        with open(self.out("Log.tsv"),"w") as file:
            file.write("Step\tMemory_MB\tTime\n")

    def __getattr__(self, name):
        return self[name]
    
    def out(self, x):
        return os.path.join(self.output_directory, x)
    
    def check_settings(self):
        # TODO - Check whether settings are of the right type 
        assert isinstance(self.seed, int)
        assert isinstance(self.classification['train_ancestry'], str)
        assert isinstance(self.classification['infer_ancestry'], str)

    def log(self, text):
        with open(self.out("Log.tsv"),"a") as file:
            file.write(
                text + 
                "\t" + 
                str(psutil.Process(os.getpid()).memory_info().rss*1.0/10**6) + 
                "\t" + 
                time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + 
                "\n")

class RScriptRunner():
    # Path to executible Rscript
    def __init__(self, r_path):
        self.r_execute = r_path
    
    def _make_executable(self, script_path):

        # Check if the script is not already executable
        if not os.access(script_path, os.X_OK):
            try:
                # Make the R script executable
                subprocess.run(['chmod', '+x', script_path], check=True)
                print(f'Permissions updated: {script_path} is now executable.')
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f'Failed to set executable permission for {script_path}: {e}')

    def run_script(self, script_path, args=None):

         # Ensure the script is executable
        self._make_executable(script_path)

        if args is not None and not isinstance(args, list):
            raise TypeError("Arguments must be provided as a list.")
        
        # Add the Rscript, script_path and arguments to the run command
        command = [self.r_execute, script_path]
        if args:
            command.extend(args)

        # Run subprocess
        subprocess.run(command)
    
    def check_rversion(self):
        # Check R version
        subprocess.run([self.r_execute, "--version"])


# Functions to normalize data
def normalize_log(X, e=0.01):
    """
    Natural log transformation.
    X: Vector or Matrix to be transformed.
    e: threshold to correct for infinity values.
    """
    X = np.log(X + e)
    return X

def normalize_minmax(X):
    """
    Columnwise min-max transforamtion. 
    X: Vector or Matrix to be transformed.
    """
    return (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))

def normalize(func, X):
    """
    Takes a normalization function to nomralize data.
    func: Normalization function.
    X: Vector or Matrix to be transformed.
    """
    normalized_data = func(X)
    return normalized_data

# Function to create integer vector from categotical column
def encode_y(data, dictionary):
    """
    Creates vector with class labels; additionally returns mapping dictionary.
    data: The data to encode.
    dictionary: Dictionary with mappings.
    """
    # Map categorical values to integers
    y = data.map(dictionary)
    y = np.array(y)
    return y


def stratified_subset(data, proportion, output_column, seed): 
    """
    Makes a subset from the training data (usually Europeans) mimicking a given frequency of classes.
    data: Data in Anndata format.
    Proportion: Frequences of classes.
    output_column: On which label to apply.
    """
    # Check for correct instances
    assert isinstance(proportion, dict), f'Proportion needs to be a dictionary'
    
    def get_sample(df, freq, seed):
        sample_size = freq[df[output_column].iloc[0]]
        return df.sample(sample_size, replace=False, random_state=seed)

    idx = data.obs.groupby(output_column, group_keys=False, observed=False).apply(lambda x: get_sample(x,proportion, seed)).index
    # Check that index is unique
    assert len(idx) == len(idx.unique())

    train_data = data[~data.obs_names.isin(idx)]
    test_data = data[data.obs_names.isin(idx)]
    return train_data, test_data

def stratified_sample(df, category_col, proportions, n):
    '''
    Take a stratified sample from a dataframe with a given number of samples and a proportion.
    df: Dataframe from which is sampled.
    category_col: Column with the classes.
    proportions: Dictionary with class proportions.
    n: Number of observations per sample.
    '''
    total_samples = df.shape[0]
    # Adjust for edge case where the total sample size is larger than available rows
    if n > total_samples:
        raise ValueError(f'Requested sample size (n={n}) exceeds the available data (n={total_samples}).')

    # Check for correct instances
    assert isinstance(df, pd.DataFrame), f'Data object needs to be pandas dataframe.'
    assert isinstance(proportion, dict), f'Proportion needs to be a dictionary'

    # Ensure proportions sum to 1 
    if sum(proportions.values()) != 1:
        raise ValueError("Proportions must sum to 1.")
    
    samples = []
    # Loop through each class 
    for category, proportion in proportions.items():
        requested_n = int(n * proportion)
        available_n = df[df[category_col] == category].shape[0]
        print(f'requested {requested_n}; available {available_n}')
        # Error occures when sample_size is bigger than population
        if requested_n > available_n:
            requested_n = available_n
        # Sample each class without replacement
        category_sample = df[df[category_col] == category].sample(n=requested_n, replace=False)
        # Append sample to list of samples
        samples.append(category_sample)
    # Concat samples and return as df
    return pd.concat(samples)

def stratified_folds(data, category_col, proportions, k):
    '''
    Trys to create equal folds but is mostly not possible.
    df: Dataframe from which is sampled.
    category_col: Column with the classes.
    proportions: Dictionary with class proportions.
    k: Number of folds.
    '''
    # Number of oberservation per fold
    n_observations = data.shape[0] // k

    folds = []
    i = 0
    while i < k:
        print(f'Fold {i}')
        # Generate fold with n_observation
        fold = stratified_sample(df=data, category_col=category_col, proportions=proportions, n=n_observations)
        # Add to list of fold
        folds.append(fold)
        # Remove sampled observations from the data
        data = data[~data.index.isin(fold.index)]
        i += 1
    # Returns folds and the data that was not used while sampling
    return folds, data
