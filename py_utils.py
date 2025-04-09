# Standard libraries
import numpy as np
import pandas as pd
import yaml
import re
import math

# Os library
import sys
import uuid
import time
from pathlib import Path
import psutil
import os
import subprocess

import pickle
import importlib

# Machine learning
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

# Data library
import anndata as ad

# Object that contains settings and methods
class Setup(dict):
    def __init__(self, config_file):
        
        # Open the default config file
        try:
            default_config = yaml.safe_load(Path("default_settings.yml").read_text())
        except yaml.YAMLError as exc:
            print(exc)

        # Open custom config file
        try:
            custom_config = yaml.safe_load(Path(config_file).read_text())
        except yaml.YAMLError as exc:
            print(exc)
        
        # Combine custom and default settings 
        final_config = {**default_config, **custom_config} 

        # Add date
        final_config["date"] = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())

        # Add id if not given
        if "id" not in final_config or not final_config["id"]:
            final_config["id"] = uuid.uuid4().hex.upper()[:10]
        
        # Print statement to inform user about their analysis
        print(f"New analysis with id: {final_config["id"]}; created: {final_config["date"]}")

        # Check required settings
        required_settings = [
            # Classification
            "output_column",
            "comparison",
            "train_ancestry",
            "infer_ancestry",
            "ancestry_column",
            # Input
            "seed",
            "data_path", 
            "data_type",
            "tech",
            # Output
            "output_directory", 
            ]
        for i in required_settings:
            assert i in final_config, f"No {i} defined but required!"  
        
        # Check that covariate has to be a list
        assert "covariate" not in final_config or isinstance(final_config["covariate"], list), \
            "If covariate exists, it must be a list."

        # Check how many classes
        if len(final_config["comparison"]) > 2:
            final_config.update({"multiclass": True})
            print(f"Multiclass problem comparing: {final_config["comparison"]}.")

        elif len(final_config["comparison"]) == 2:
            final_config.update({"multiclass": False})
            print(f"Binaryclass problem comparing: {final_config["comparison"]}.")
        
        else:
            print(f"Cant do a classification with one assigned class in comparison: {final_config["classification"]["comparison"]}.")

        # Init
        self.final_config = final_config
        super().__init__(final_config)

        # Check all the settings
        self.check_settings()

    def __getattr__(self, name):
        return self[name]

    def make_output_directory(self):
        # Create output directory
        Path(self.final_config["output_directory"]).mkdir(parents=True, exist_ok=self.overwrite)
        save_location = os.path.join(os.getcwd(), self.final_config["output_directory"])
        # Print tatement to inform about save location
        print(f"Output will be saved to {save_location}")

        # Create a log file
        with open(self.out("Log.tsv"),"w") as file:
            file.write("Step\tMemory_MB\tTime\n")
        
    def out(self, x):
        return os.path.join(self.output_directory, x)
    
    def check_settings(self):
        # TODO - Check whether settings are of the right type 
        assert isinstance(self.seed, int)

        # Check data
        assert isinstance(self.data_path, str),     f"{self.data_path} needs to be a string."
        assert os.path.exists(self.data_path),      f"{self.data_path} does not exist."
        assert self.data_path.endswith(".h5ad"),    f"Only support data files in h5ad format."
        
        # Type
        assert self.tech in ["methylation", "transcriptomics", "proteomics"], \
            f"Invalid data type: {self.data_type}."


        # Classification
        assert isinstance(self.ancestry_column, str)
        assert isinstance(self.train_ancestry, str)
        assert isinstance(self.infer_ancestry, str)

    
        # Machine learning
        assert isinstance(self.nfolds, int)
        assert self.nfolds >= 2, \
            f"Cross-validation requires at least nfolds: 2 folds, got nfolds: {self.nfolds}."
        
        # Grid search
        # assert all(isinstance(x, float) for x in self.grid_search['l1_ratio']), \
        #     f"Grid search l1_ratio needs to be floating points between 0.0 and 1.0."

    def log(self, text):
        with open(self.out("Log.tsv"),"a") as file:
            file.write(
                text + 
                "\t" + 
                str(psutil.Process(os.getpid()).memory_info().rss*1.0/10**6) + 
                "\t" + 
                time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + 
                "\n")

    def return_settings(self):
        return self.final_config

    def add_settings(self, new_settings: dict, overwrite=False):
        """
        Add new key-value pairs to the settings.

        Parameters:
        -----------
        new_settings : dict
            A dictionary containing the new settings to add.
        overwrite : bool, optional (default=False)
            If True, existing settings with the same key will be replaced.
            If False, a warning will be printed if a key already exists.
        """
        for key, value in new_settings.items():
            if key in self and not overwrite:
                print(f"Setting '{key}' already exists. Use `overwrite=True` to modify it.")
            else:
                self.final_config[key] = value  # Update internal dictionary
                self[key] = value  # Ensure settings are accessible as attribute    
               
class ErrorHandler():
    def __init__(self, log_file):
        """
        Initializes the error handler.

        :param log_file: Path to the log file.
        """
        self.log_file = log_file

    def handle_error(self, error):
        """
        Handles an exception by logging it, printing it, and terminating the script.

        :param error: The exception instance to handle.
        """
        error_message = f"Error: {str(error)}"

        # Log the error
        self._log_to_file(error_message)

        # If running in an interactive session, do not terminate the session
        if not self._is_interactive():
            raise error  # Raise the error to maintain stack trace
            sys.exit(1)  # Terminate the script with a non-zero exit code

        raise error

    def _is_interactive(self):
        """
        Determines if the script is being run in an interactive environment.

        :return: True if running interactively, False otherwise.
        """
        return hasattr(sys, 'ps1') or 'IPython' in sys.modules

    def _log_to_file(self, message):
        """
        Logs a message to the specified log file.

        :param message: The message to log.
        """
        try:
            with open(self.log_file, "a") as file:
                file.write(message + "\n")
        except Exception as log_error:
            print(f"Failed to write to log file: {str(log_error)}", file=sys.stderr)

class DataValidator():
    """
    Checks if data and settings are compatible.
        
    :param data: AnnData object containing the data.
    :param setup: Setup object containing settings.
    :param error_handler: An instance of ErrorHandler to manage errors.
    """
    def __init__(self, data: ad.AnnData, setup: Setup, error_handler: ErrorHandler):
        self.data = data
        self.setup = setup
        self.error_handler = error_handler

    def validate_features(self):
        """
        Check that features are unique.
        """
        try:
            assert self.data.var.index.shape[0] == len(set(self.data.var.index)), \
                'Features not unique.'
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def validate_output_column(self):
        """
        Check that output column is in the data.
        """
        try:
            output_column = self.setup.output_column
            assert output_column in self.data.obs.columns, \
                f"Output column '{output_column}' not in the data."
            
            has_na = self.data.obs[output_column].isnull().any()
            assert not has_na, \
                f"The output column '{output_column}' contains NA values. Please handle missing data before proceeding."
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def validate_class_labels(self):
        """
        Check that class labels of comparisons are present in the output column of the data.
        """
        try:
            required_labels = self.setup.comparison
            output_column = self.setup.output_column
            
            for label in required_labels:
                assert self.data.obs[output_column].str.contains(label).any(), \
                    f"No '{label}' in output column: '{output_column}', choose different comparison in settings."
        except AssertionError as e:
            self.error_handler.handle_error(e)
            
    def validate_ancestry_column(self):
        """
        Check if the ancestry column is present in the data.
        """
        try:
            ancestry_column = self.setup.ancestry_column
            assert ancestry_column in self.data.obs.columns, \
                f"Ancestry column '{ancestry_column}' not in the data."
            
            has_na = self.data.obs[ancestry_column].isnull().any()
            assert not has_na, \
                f"The ancestry column '{ancestry_column}' contains NA values. Please handle missing data before proceeding."
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def validate_ancestry(self):
        """
        Check if ancestries are present in ancestry column.
        """
        try:
            ancestry_column = self.setup.ancestry_column
            train_ancestry = self.setup.train_ancestry
            inf_ancestry = self.setup.infer_ancestry

            assert self.data.obs[ancestry_column].str.contains(train_ancestry).any(), \
                f"No '{train_ancestry}' in ancestry_column: '{ancestry_column}'."
            assert self.data.obs[ancestry_column].str.contains(inf_ancestry).any(), \
                f"No '{inf_ancestry}' in ancestry_column: '{ancestry_column}'."
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def validate_covariate_column(self):
        """
        Check if covariate column(s) are present in the data, if specified.
        Check that there are no NA values in the covariate column(s).
        """
        try:
            covariates = self.setup.get("covariate")
            if covariates:  
                if not isinstance(covariates, list):
                    covariates = [covariates]  
                
                for covariate in covariates:
                    assert covariate in self.data.obs.columns, \
                        f"No '{covariate}' column in the data, can't be used as covariate."
                    
                    has_na = self.data.obs[covariate].isnull().any()
                    assert not has_na, \
                        f"The covariate column '{covariate}' contains NA values. Please handle missing data before proceeding."
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def validate_na_counts(self):
        """
        Check if there are NA values in the count matrix (.X) and provide detailed information, 
        including the proportion of NA values in the error message.
        """
        try:
            if hasattr(self.data.X, "toarray"):  
                dense_matrix = self.data.X.toarray()  # Convert sparse matrix to dense if needed
            else:
                dense_matrix = self.data.X  # Already a dense matrix

            # Calculate the total number of elements in the matrix
            total_elements = dense_matrix.size

            # Count how many NA (NaN) values are in the matrix
            na_elements = np.sum(np.isnan(dense_matrix))

            # Calculate the proportion of NA values
            proportion_na = na_elements / total_elements

            # Identify rows and columns with NA values
            rows_with_all_na = np.where(np.isnan(dense_matrix).all(axis=1))[0]
            obs_with_all_na = self.data.obs.index[rows_with_all_na].tolist()

            cols_with_any_na = np.where(np.isnan(dense_matrix).any(axis=0))[0]
            var_with_any_na = self.data.var.index[cols_with_any_na].tolist()

            # Raise error if there are NA values, including the proportion of NA values in the message
            if len(obs_with_all_na) > 0 or len(var_with_any_na) > 0:
                raise ValueError(
                    f"NA values found in the count matrix. Proportion of NA values: {proportion_na:.4f}. "
                    f"There are {len(obs_with_all_na)} observations with all NA values and "
                    f"{len(var_with_any_na)} columns with any NA values.\n"
                    f"Please check the data preprocessing steps."
                )

        except ValueError as e:
            self.error_handler.handle_error(e)
    
    def validate_negative_counts(self):
        """
        Check if there are any negative values in the count matrix (.X).
        Raise an error with the proportion of negative values if present.
        """
        try:
            if hasattr(self.data.X, "toarray"):  
                dense_matrix = self.data.X.toarray()  # Convert sparse matrix to dense if needed
            else:
                dense_matrix = self.data.X  # Already a dense matrix

            # Calculate the total number of elements in the matrix
            total_elements = dense_matrix.size

            # Count how many negative values are in the matrix
            negative_elements = np.sum(dense_matrix < 0)

            # Calculate the proportion of negative values
            proportion_negative = negative_elements / total_elements

            # Raise error if there are any negative values and include the proportion in the error message
            if negative_elements > 0:
                raise ValueError(f"Negative values found in the count matrix. " 
                                f"Proportion of negative values: {proportion_negative:.4f}.\n"
                                f"Please check the data preprocessing steps.")
        
        except ValueError as e:
            self.error_handler.handle_error(e)


    def check_min_samples_per_class(self, data, column, min_samples, data_name):
        """
        Check if classes in a specified column have at least a minimum number of samples.
        """
        try:
            counts = data[column].value_counts()
            failing_classes = counts[counts < min_samples]
            assert failing_classes.empty, (
                f"Minimum sample size of {min_samples} per class is required. "
                f"The following classes in '{data_name}' fail to meet the requirement:\n{counts}"
            )
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def check_data_leakage(self, train_idx, test_idx):
        """
        Check for data leakage between training and test sets.
        """
        try:
            overlapping_indices = test_idx[test_idx.isin(train_idx)]
            assert overlapping_indices.empty, (
                "Data leakage detected! Some observations in the test set also occur in the training set."
            )
        except AssertionError as e:
            self.error_handler.handle_error(e)

    def data_settings_compatibility(self):
        """
        Run all validation checks on the data.
        """
        try:
            self.validate_features()
            self.validate_output_column()
            self.validate_class_labels()
            self.validate_ancestry_column()
            self.validate_ancestry()
            self.validate_covariate_column()
        except AssertionError as e:
            self.error_handler.handle_error(e)

class ScriptRunner():
    """
    Initialize ScriptRunner with paths to R and Python executables.
        
    :param r_path: Path to the R executable (default is 'Rscript')
    :param py_path: Path to the Python executable (default is 'python')
    """

    # Path to executible Rscript
    def __init__(self, r_path, py_path):
        self.r_execute = r_path
        self.py_execute = py_path
    
    def _make_executable(self, script_path):
        '''
        Ensures that the script is executable.
        '''

        # Check if the script is not already executable
        if not os.access(script_path, os.X_OK):
            try:
                # Make the R script executable
                subprocess.run(['chmod', '+x', script_path], check=True)
                print(f'Permissions updated: {script_path} is now executable.')
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f'Failed to set executable permission for {script_path}: {e}')

    def run_script(self, script_path, args=None):

        # Determine script type
        if script_path.endswith('.R'):
            executable = self.r_execute
        elif script_path.endswith('.py'):
            executable = self.py_execute
        else:
            raise ValueError('Script type not supported. Please provide an .R or .py file.')

        # Ensure the script is executable
        self._make_executable(script_path)

        if args is not None and not isinstance(args, list):
            raise TypeError('Arguments must be provided as a list.')
        
        # Add the Rscript, script_path and arguments to the run command
        command = [executable, script_path]
        if args:
            command.extend(args)

        # Run subprocess
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Script {script_path} failed with exit code {e.returncode}")
    
    def check_rversion(self):
        # Check R version
        subprocess.run([self.r_execute, "--version"])

# Normalization
def normalize_log(X, e=1):
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

def beta_to_mvalue(X, epsilon=0.00001):
    """
    Convert beta values to M-values using the formula:
    M = log2(X / (1 - X)), with clipping to avoid extreme values.

    Parameters:
    X (array-like or float): X values (0 < X < 1)
    epsilon (float): Small value to prevent log(0) or division by zero (0 <= epsilon <= 0.5)

    Returns:
    numpy.ndarray or float: M-values
    """
    # Validate input types
    if not isinstance(X, (list, np.ndarray, float, int)):
        raise ValueError("Invalid value for X: Expected numeric values or an array.")

    if not isinstance(epsilon, (float, int)) or not (0 <= epsilon <= 0.5):
        raise ValueError("Invalid value for epsilon; expected 0 <= epsilon <= 0.5")

    # Convert to NumPy array for vectorized operations
    X = np.array(X, dtype=np.float64)

    # Clip X values to avoid extreme cases
    X = np.clip(X, epsilon, 1 - epsilon)

    # Convert X values to M-values
    m_values = np.log2(X / (1 - X))

    return m_values

def raw(X):
    return(X)

# Dictionary to look up methods to normalize data
ml_normalization_methods = {
    "transcriptomics": {
        "normalize_log": normalize_log,
        "normalize_minmax": normalize_minmax,
        "raw": raw
    },
    "methylation": {
        "beta_to_mvalue": beta_to_mvalue,
        "normalize_log": normalize_log,
        "normalize_minmax": normalize_minmax,
        "raw": raw
    },
    "proteomics": {
        "raw": raw
    }
}

# Filtering features
def calculate_tmm_norm_factors(data):
    """
    Calculates TMM normalization factors for RNA-Seq data.
    
    Parameters:
    - data: A numpy array with samples as rows and genes as columns (samples x genes).
    
    Returns:
    - A numpy array containing the normalization factors for each sample.
    """
    # Step 1: Compute the geometric mean of the counts for each gene (column)
    geometric_mean = np.exp(np.mean(np.log(data + 1), axis=0))  # Axis 0: Across samples

    # Step 2: Compute M-values (log-ratio) between each sample and the geometric mean
    M_values = data / geometric_mean  # Element-wise division by geometric mean
    M_values = np.log2(M_values + 1)  # Add pseudocount for stability
    
    # Step 3: Trim the extreme M-values to reduce the influence of highly differentially expressed genes
    trim_percent = 0.05
    def trim_m_values(x):
        lower = np.quantile(x, trim_percent)
        upper = np.quantile(x, 1 - trim_percent)
        x = np.clip(x, lower, upper)  # Clip values to the trim range
        return x
    
    # Apply trimming to each sample (row-wise)
    M_values_trimmed = np.apply_along_axis(trim_m_values, axis=1, arr=M_values)

    # Step 4: Calculate the normalization factors
    norm_factors = np.median(M_values_trimmed, axis=1)  # Median across genes for each sample

    # Normalize by the median normalization factor
    norm_factors = norm_factors / np.median(norm_factors)
    
    return norm_factors

# CPM
def cpm(data, norm_factors=None, log=False):
    """
    Converts raw counts to CPM (Counts Per Million), optionally using normalization factors.
    Optionally, log-transform the CPM values.

    Parameters:
    - data: A numpy array with rows as samples and columns as genes (samples x genes).
    - norm_factors: A numpy array with normalization factors for each sample (length = n_samples). Default is None.
    - log: Boolean value, whether to log-transform the CPM values. Default is False.

    Returns:
    - A numpy array with CPM or log-transformed CPM values (samples x genes).
    """
    data = np.asarray(data)  # Ensure data is a NumPy array

    # Check that norm_factors matches the number of samples
    if norm_factors is not None:
        norm_factors = np.asarray(norm_factors)
        if norm_factors.shape[0] != data.shape[0]:
            raise ValueError("The number of samples must match the length of the normalization factors.")

        # Apply normalization factors (element-wise division across rows)
        data = data / norm_factors[:, np.newaxis]  # Broadcast over columns

    # Calculate total reads per sample (sum of counts for each sample)
    total_reads_per_sample = np.sum(data, axis=1) + 1e-6  # Prevent division by zero

    # Convert raw counts to CPM
    cpm_data = (data / total_reads_per_sample[:, np.newaxis]) * 1e6  # Broadcast over columns

    # Apply log transformation if requested
    if log:
        cpm_data = np.log2(cpm_data + 1)  # Add pseudocount for stability

    return cpm_data

# By signal
def signal_by_percentile(data, percentile):
    """
    Calculates the signal threshold for genes based on the specified percentile of their total counts.
    
    Parameters:
    - data: A pandas DataFrame with rows as samples and columns as genes (samples x genes).
    - percentile: The percentile to calculate the signal threshold (e.g., 25 for 25th percentile).
    
    Returns:
    - A numeric value representing the calculated signal threshold for genes based on the specified percentile.
    """
    # Calculate total signal for each gene across all samples (sum along columns)
    gene_signal = data.sum(axis=0)
    
    # Calculate the threshold using the specified percentile
    signal_threshold = np.percentile(gene_signal, percentile)
    
    return signal_threshold

def filter_by_signal(data, min_signal=10, max_signal=None, min_samples_ratio=0.5):
    """
    Filters genes based on minimum and maximum signal and presence in at least a given fraction of samples.
    
    Parameters:
    - data: A pandas DataFrame with rows as samples and columns as genes (samples x genes).
    - min_signal: The minimum total signal required for a gene (default: 10).
    - max_signal: The maximum total signal for a gene (optional).
    - min_samples_ratio: The minimum fraction of samples a gene must be expressed in (default: 50%).
    
    Returns:
    - A boolean pandas Series indicating which genes were retained.
    """
    # Calculate the number of samples
    num_samples = data.shape[0]
    
    # Calculate the minimum number of samples for the given ratio
    min_samples = int(np.floor(min_samples_ratio * num_samples))
    
    # Calculate total signal for each gene across all samples (sum along columns)
    gene_signal = data.sum(axis=0)
    
    # Calculate how many samples express each gene (non-zero signal counts)
    gene_samples = (data > 0).sum(axis=0)
    
    # Filter genes based on min_signal, max_signal, and min_samples_ratio
    if max_signal is None:
        # If max_signal is not provided, only filter by min_signal and sample count
        filtered_genes = (gene_signal >= min_signal) & (gene_samples >= min_samples)
    else:
        # If max_signal is provided, filter by both min_signal and max_signal
        filtered_genes = (gene_signal >= min_signal) & (gene_signal <= max_signal) & (gene_samples >= min_samples)
    
    return filtered_genes

# By variance
def variance_by_percentile(data, percentile=25):
    """
    Calculates the variance threshold for genes based on the specified percentile.
    
    Parameters:
    - data: A pandas DataFrame with rows as samples and columns as genes (samples x genes).
    - percentile: The percentile to calculate the variance threshold (default: 25).
    
    Returns:
    - A numeric value representing the calculated variance threshold based on the specified percentile.
    """
    # Compute variance for each gene across samples (variance along columns)
    gene_variances = data.var(axis=0)
    
    # Determine the variance threshold at the given percentile
    variance_threshold = np.percentile(gene_variances, percentile)
    
    return variance_threshold

def filter_by_variance(data, var_threshold=0.01, min_samples_ratio=0.5):
    """
    Filters methylation genes based on variance and presence in at least a given fraction of samples.
    
    Parameters:
    - data: A pandas DataFrame with rows as samples and columns as genes (samples x genes).
    - var_threshold: The minimum variance threshold to retain a gene (default: 0.01).
    - min_samples_ratio: The minimum fraction of samples where the gene must be methylated above 0 (default: 50%).
    
    Returns:
    - A boolean pandas Series indicating which genes were retained.
    """
    # Calculate the number of samples
    num_samples = data.shape[0]
    
    # Calculate the minimum number of samples for the given ratio
    min_samples = int(np.floor(min_samples_ratio * num_samples))
    
    # Calculate variance for each gene across samples (variance along columns)
    gene_variances = data.var(axis=0)
    
    # Count how many samples have a methylation value > 0 for each gene (non-zero beta value)
    methylated_samples = (data > 0).sum(axis=0)
    
    # Filter genes based on variance threshold and the number of methylated samples
    filtered_genes = (gene_variances > var_threshold) & (methylated_samples >= min_samples)
    
    return filtered_genes



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

def classify_covariates(df, covariates):
    """
    Classifies the given covariates into continuous and discrete groups.
    
    Parameters:
    - df: The DataFrame containing the data.
    - covariates: List of covariates (column names) to classify.
    
    Returns:
    - A dictionary with two keys: 'continuous' and 'discrete', each containing a list of covariates.
    """
    classified_covariates = {'continuous': [], 'discrete': []}
    
    for covariate in covariates:
        if covariate not in df.columns:
            continue  # Skip covariates not found in the DataFrame
        
        col_data = df[covariate]
        unique_values = col_data.nunique()
        
        # Check for categorical dtype
        if pd.api.types.is_categorical_dtype(col_data):
            classified_covariates['discrete'].append(covariate)
        elif np.issubdtype(col_data.dtype, np.number):
            # Heuristic: Consider numerical columns with few unique values as discrete
            if unique_values <= 10:  # You can adjust this threshold
                classified_covariates['discrete'].append(covariate)
            else:
                classified_covariates['continuous'].append(covariate)
        else:
            # Non-numerical columns (e.g., strings) are discrete
            classified_covariates['discrete'].append(covariate)
    
    return classified_covariates

def stratified_subset(data, freq_dict, group_columns, seed):
    """
    Sample from the DataFrame based on the provided sample_sizes dictionary.
    
    Parameters:
    - df: DataFrame with the data to sample from.
    - sample_sizes: Dictionary with the sample sizes for each unique combination of the group_columns.
    - group_columns: List of column names to group by.
    
    Returns:
    - A DataFrame containing the sampled rows.
    """
    
    # Function to sample based on the provided dictionary
    def sample_groups(group, sample_sizes, group_columns, seed):
            # Create a key based on the group columns
            if len(group_columns) == 1:
                group_key = group[group_columns[0]].iloc[0]  # Scalar for single-column group
            else:
                group_key = tuple(group[col].iloc[0] for col in group_columns)  # Tuple for multi-column group
            
            # Get the sample size for this group
            sample_size = sample_sizes.get(group_key, 0)
            
            if sample_size > 0:
                return group.sample(n=sample_size, replace=False, random_state=seed)  # Sample without replacement
            else:
                return pd.DataFrame() 

    # Apply the sampling function to each group based on the dictionary and group_columns
    idx = data.obs.groupby(group_columns, group_keys=False, observed=False).apply(sample_groups, sample_sizes=freq_dict, group_columns=group_columns, seed=seed).index
    # Check that index is unique
    assert len(idx) == len(idx.unique())

    train_data = data[~data.obs_names.isin(idx)]
    test_data = data[data.obs_names.isin(idx)]
    return train_data, test_data


# def stratified_subset(data, proportion, output_column, seed): 
#     """
#     Makes a subset from the training data (usually Europeans) mimicking a given frequency of classes.
#     data: Data in Anndata format.
#     Proportion: Frequences of classes.
#     output_column: On which label to apply.
#     """
#     # Check for correct instances
#     assert isinstance(proportion, dict), f'Proportion needs to be a dictionary'
    
#     def get_sample(df, freq, seed):
#         sample_size = freq[df[output_column].iloc[0]]
#         return df.sample(sample_size, replace=False, random_state=seed)

#     idx = data.obs.groupby(output_column, group_keys=False, observed=False).apply(lambda x: get_sample(x,proportion, seed)).index
#     # Check that index is unique
#     assert len(idx) == len(idx.unique())

#     train_data = data[~data.obs_names.isin(idx)]
#     test_data = data[data.obs_names.isin(idx)]
#     return train_data, test_data

def train_algorithm(algo_name, setup, train_X, train_y, available_jobs, prop=None):
    """
    Train the model using GridSearchCV and save the model and hyperparameters.
    All parameters need to be pickable.

    Parameters:
        algo_name (str): Name of the algorithm to be trained.
        setup (dict): Dictionary with settings.
        train_X (np.array): Training feature matrix.
        train_y (np.array): Training target vector.
        available_jobs (int): Number of CPUs to use.

    Returns:
        best_model (sklearn model): The best-trained model from GridSearchCV.
    """
    try:

        # Dynamically load the appropriate scikit-learn module based on algorithm name
        module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
        algo_class = getattr(importlib.import_module(module_name), algo_name)

        # Create an instance of the algorithm with the provided random seed
        algo_instance = algo_class(random_state=setup["seed"])

        # Define cross-validation strategy
        cv_hyp = StratifiedKFold(n_splits=setup["nfolds"], shuffle=True, random_state=setup["seed"])

        # Load hyperparameter search space for the algorithm
        search_space = setup["grid_search"][algo_name]

        # Perform hyperparameter tuning using GridSearchCV
        grid_search_jobs = max(1, math.floor(available_jobs))
        grid_search = GridSearchCV(
            estimator=algo_instance,
            param_grid=search_space,
            cv=cv_hyp,
            scoring="f1_weighted",
            refit=True,
            n_jobs=grid_search_jobs,  # Parallel jobs for GridSearchCV
        )
        grid_search.fit(train_X, train_y)  # Train the model on the training data
        best_model = grid_search.best_estimator_  # Extract the best model from the grid search

        # Save the model
        if setup["save_model"]:
            if prop is not None:
                with open(os.path.join(setup["output_directory"], f"{algo_name}_{prop}_{setup["id"]}.pkl"), "wb") as f:
                    pickle.dump(best_model, f, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                  with open(os.path.join(setup["output_directory"], f"{algo_name}_{setup["id"]}.pkl"), "wb") as f:
                    pickle.dump(best_model, f, protocol=pickle.HIGHEST_PROTOCOL)


        # Training complete
        print(f"{algo_name} training done.")
        return best_model

    except Exception as e:
        # Handle and log any exceptions that occur during training
        print(f"Error occurred during training for {algo_name}: {e}")
        return None
    
def evaluate_algorithm_cross_ancestry(
        algo_name, 
        best_model, 
        setup, 
        test_X, 
        test_y, 
        inf_X, 
        inf_y, 
        feature_names,
        encoder,
        prop=None
        ):
    """
    Evaluate the trained model on test and inference datasets, saving results and feature importances.

    Parameters:
        algo_name (str): Name of the algorithm being evaluated.
        setup (object): Custom setup object containing configuration like grid search, output paths, and seeds.
        best_model (sklearn model): The trained model to be evaluated.
        test_X (np.array): Test feature matrix.
        test_y (np.array): Test target vector.
        inf_X (np.array): Inference feature matrix for a separate subset.
        inf_y (np.array): Inference target vector.
        feature_names (list): List of feature names.
        encoder (dict): Dictionary of class names.

    Returns:
        Probabilities, Interpretations and Hyperparameters
    """
    try:
        # Log the start of testing
        setup.log("Testing")

        # Generate predictions on the test set and save probabilities
        test_y_hat = pd.DataFrame(best_model.predict_proba(test_X))
        test_y_hat = test_y_hat.rename(columns=encoder)
        test_y_hat["y"] = test_y
        # Additional information
        test_y_hat["Status"] = "Test"
        test_y_hat["Prediction"] = "Subset"

        # Log the start of inference
        setup.log("Infering")

        # Generate predictions on the inference set and save probabilities
        inf_y_hat = pd.DataFrame(best_model.predict_proba(inf_X))
        inf_y_hat = inf_y_hat.rename(columns=encoder)
        inf_y_hat["y"] = inf_y
        # Additional information
        inf_y_hat["Status"] = "Inference"
        inf_y_hat["Prediction"] = "Ancestry"
        
        # Combine propabilities (test inference)
        y_hat = pd.concat([test_y_hat, inf_y_hat])
        y_hat["Seed"] = setup.seed
        y_hat["Algorithm"] = algo_name

        # Inference complete
        print(f"{algo_name} validation done.")

        # Log the start of model interpretation
        setup.log("Model interpretations")
        # Extract feature importance
        feature_importance = extract_feature_importance(best_model, feature_names)
        feature_importance["Seed"] = setup.seed
        feature_importance["Algorithm"] = algo_name

        # Log the start of hyperparameters
        hyperparameters = best_model.get_params()
        hyperparameters = pd.DataFrame([hyperparameters])
        hyperparameters["Seed"] = setup.seed
        hyperparameters["Algorithm"] = algo_name

        return y_hat, feature_importance, hyperparameters

    except Exception as e:
        # Handle and log any exceptions that occur during evaluation
        print(f"Error occurred during evaluation for {algo_name}: {e}")

def evaluate_algorithm_eur_subsetting(
        algo_name, 
        best_model, 
        subset_X, 
        subset_y
        ):
    """
    Evaluate the trained model on test subsets and return probabilities and metrics.
    """
    try:
        propability_list = []
    
        # Iterate over all subsets (X, y)
        for X, y in zip(subset_X, subset_y):
            # Get predicted probabilities
            test_y_hat = pd.DataFrame(best_model.predict_proba(X))
            test_y_hat["y"] = y
            test_y_hat["Algorithm"] = algo_name
            test_y_hat["n_test_ancestry"] = y.shape[0]
            
            # Collect probabilities
            propability_list.append(test_y_hat)

        # Combine probabilities
        propabilities = pd.concat(propability_list, ignore_index=True)

        # Finished validation
        print(f"{algo_name} validation done.")

        return propabilities
    
    except Exception as e:
        print(f"Error occurred during evaluation for {algo_name}: {e}")


def evaluate_algorithm_robustness(
        algo_name, 
        best_model, 
        setup, 
        test_X, 
        test_y, 
        inf_X, 
        inf_y, 
        feature_names, 
        prop,
        encoder 
        ):
    """
    Evaluate the trained model, save probabilities, calculate metrics, and return results.
    """
    try:
        # Log the start of testing
        setup.log("Testing")

        # # Generate predictions on the test set and save probabilities
        test_y_hat = pd.DataFrame(best_model.predict_proba(test_X))
        test_y_hat = test_y_hat.rename(columns=encoder)
        test_y_hat["y"] = test_y
        # Additional information
        test_y_hat["Algorithm"] = algo_name
        test_y_hat["Status"] = "Test"
        test_y_hat["Prediction"] = "Subset"
        # test_y_hat.to_csv(setup.out(f"Probabilities_{prop}_{algo_name}_test.csv"), index=False)

        # Log the start of inference
        setup.log("Infering")

        # Generate predictions on the inference set and save probabilities
        inf_y_hat = pd.DataFrame(best_model.predict_proba(inf_X))
        inf_y_hat = inf_y_hat.rename(columns=encoder)
        inf_y_hat["y"] = inf_y
        # Additional information
        inf_y_hat["Algorithm"] = algo_name
        inf_y_hat["Status"] = "Inference"
        inf_y_hat["Prediction"] = "Ancestry"
        inf_y_hat["y"] = inf_y
        # inf_y_hat.to_csv(setup.out(f"Probabilities_{prop}_{algo_name}_inf.csv"), index=False)

        # Combine propabilities
        y_hat = pd.concat([test_y_hat, inf_y_hat])
        y_hat["Seed"] = setup.seed
        # Save
        # y_hat.to_csv(setup.out(f"Probabilities_{algo_name}.csv"), index=False)

        print(f"{algo_name} validation done.")
        
        # Feature importance extraction
        setup.log("Model interpretations")
        feature_importance = extract_feature_importance(best_model, feature_names)
        feature_importance.to_csv(setup.out(f"Feature_importance_{prop}_{algo_name}.csv"), index=False)

        # Calculate metrics
        test_probabilities, test_y = test_y_hat.iloc[:, 1].values, test_y_hat["y"].values
        inf_probabilities, inf_y = inf_y_hat.iloc[:, 1].values, inf_y_hat["y"].values

        # Calculate AUC scores
        test_auc_score = roc_auc_score(y_true=test_y, y_score=test_probabilities)
        inf_auc_score = roc_auc_score(y_true=inf_y, y_score=inf_probabilities)

        # Prepare dataframes for metrics
        train_data_metric = {algo_name: test_auc_score, "Status": "Test", "Prediction": "Subset"}
        inf_data_metric = {algo_name: inf_auc_score, "Status": "Inference", "Prediction": "Ancestry"}
        
        # Combine metrics
        test_metric = pd.DataFrame(train_data_metric, index=[0])
        inf_metric = pd.DataFrame(inf_data_metric, index=[0])
        metric_df = pd.concat([test_metric, inf_metric])

        return y_hat

    except Exception as e:
        print(f"Error occurred during evaluation for {algo_name}: {e}")


# def run_algorithm(algo_name, setup, train_X, train_y, test_X, test_y, inf_X, inf_y, feature_names, available_jobs):
#     """
#     Train, validate, and save results for a given algorithm in parallel.

#     Parameters:
#         algo_name (str): Name of the algorithm to be used (e.g., 'LogisticRegression').
#         setup (object): Custom setup object containing configuration like grid search, output paths, and seeds.
#         train_X (np.array): Training feature matrix.
#         train_y (np.array): Training target vector.
#         test_X (np.array): Test feature matrix.
#         test_y (np.array): Test target vector.
#         inf_X (np.array): Inference feature matrix for a separate subset.
#         inf_y (np.array): Inference target vector.
#         feature_names (list): List of feature names.
#         available_jobs (int): Number of cpus to use.

#     Returns:
#         None
#     """
#     try:
#         # Log the algorithm being processed
#         setup.log(algo_name)

#         # Dynamically load the appropriate scikit-learn module based on algorithm name
#         module_name = "sklearn.linear_model" if "LogisticRegression" in algo_name else "sklearn.ensemble"
#         algo_class = getattr(importlib.import_module(module_name), algo_name)

#         # Create an instance of the algorithm with the provided random seed
#         algo_instance = algo_class(random_state=setup.seed)

#         # Define cross-validation strategy
#         cv_hyp = StratifiedKFold(n_splits=setup.nfolds, shuffle=True, random_state=setup.seed)

#         # Load hyperparameter search space for the algorithm
#         search_space = setup.grid_search[algo_name]

#         # Log the start of training
#         setup.log("Training")

#         # Perform hyperparameter tuning using GridSearchCV
#         grid_search_jobs = max(1, math.floor(available_jobs))
#         grid_search = GridSearchCV(
#             estimator=algo_instance,
#             param_grid=search_space,
#             cv=cv_hyp,
#             scoring="f1_weighted",
#             refit=True,
#             n_jobs=grid_search_jobs,  # Parallel jobs for GridSearchCV
#         )
#         grid_search.fit(train_X, train_y)  # Train the model on the training data
#         best_m = grid_search.best_estimator_  # Extract the best model from the grid search

#         # Save the best hyperparameters to a CSV file
#         hyp = pd.DataFrame(grid_search.best_params_, index=[setup.id])
#         hyp.to_csv(setup.out(f"Hyperparameters_{algo_name}.csv"), index=False)

#         # Optionally save the trained model to a pickle file
#         if setup.save_model:
#             with open(setup.out(f"{algo_name}_{setup.id}.pkl"), "wb") as f:
#                 pickle.dump(best_m, f, protocol=pickle.HIGHEST_PROTOCOL)

#         # Training complete
#         print(f"{algo_name} training done.")

#         # Log the start of testing
#         setup.log("Testing")

#         # Generate predictions on the test set and save probabilities
#         test_y_hat = pd.DataFrame(best_m.predict_proba(test_X))
#         test_y_hat["y"] = test_y
#         test_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_test.csv"), index=False)

#         # Log the start of inference
#         setup.log("Infering")

#         # Generate predictions on the inference set and save probabilities
#         inf_y_hat = pd.DataFrame(best_m.predict_proba(inf_X))
#         inf_y_hat["y"] = inf_y
#         inf_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_inf.csv"), index=False)

#         # Inference complete
#         print(f"{algo_name} validation done.")

#         # Log the start of model interpretation
#         setup.log("Model interpretations")

#         # Extract feature importances (if supported by the model)
#         feature_importance = extract_feature_importance(best_m, feature_names)

#         # Save feature importances to a CSV file
#         feature_importance.to_csv(setup.out(f"Feature_importance_{algo_name}.csv"), index=False)

#     except Exception as e:
#         # Handle and log any exceptions that occur during the process
#         print(f"Error occurred while processing {algo_name}: {e}")



# A function to extract feature importance
def extract_feature_importance(model, feature_names):
    """
    Extract feature importance from a given model.
    
    Parameters:
        model: Trained sklearn model
        feature_names: List of feature names (columns)
    
    Returns:
        DataFrame with feature importance scores
    """
    if hasattr(model, "feature_importances_"):  # Tree-based models
        importance = model.feature_importances_
    elif hasattr(model, "coef_"):  # Linear models
        importance = model.coef_.flatten()
    else:
        raise ValueError(f"Model {type(model).__name__} does not support feature importance extraction.")
    
    # Normalize importance scores to sum to 1 (optional)
    importance = importance / np.sum(importance)
    
    # Create a DataFrame for readability
    importance_df = (
        pd.DataFrame({"Feature": feature_names, "Importance": importance})
        .sort_values(by="Importance", ascending=False)
        .reset_index(drop=True)
        )
    
    return importance_df


# EUR Subsetting / Robustness
def sample_by_size(adata, props, seed, output_column):
    """
    Sample subsets of data based on proportions while ensuring that each class has at least two samples.

    :param adata: The data to sample from (AnnData or similar).
    :param props: The list of proportions to sample.
    :param seed: Random seed for reproducibility.
    :param output_column: The column to group by for stratification.
    :param error_handler: An instance of the ErrorHandler class to manage errors.
    """
    # Set seed for reproducibility
    np.random.seed(seed)  
    
    # List to store the sampled subsets
    samples = {}
    
    # Ensure each class has at least two samples
    class_groups = adata.obs.groupby(output_column, observed=False)
    
    for prop in props:
        # Name for the subset
        sample_name = f"proportion_{str(prop).replace('.', '_')}"
        
        # Calculate the total number of samples for this proportion
        total_size = int(len(adata) * prop)
        
        # Initialize storage for indices
        selected_indices = []
        
        # First, ensure at least 2 samples per class
        for class_name, class_data in class_groups:
            if len(class_data) >= 2:
                selected_indices += class_data.sample(n=2, random_state=seed).index.tolist()
        
        # Remaining size to sample after ensuring 2 samples per class
        remaining_size = total_size - len(selected_indices)
        
        if remaining_size > 0:
            # Create a DataFrame excluding the already selected indices
            remaining_data = adata[~adata.obs.index.isin(selected_indices)]
            
            # Randomly sample remaining indices
            additional_indices = remaining_data.obs.sample(n=remaining_size, random_state=seed).index.tolist()
            selected_indices += additional_indices
        
        # Store the sampled subset
        samples[sample_name] = adata[selected_indices]
    
    return samples

# # EUR Subsetting / Robustness
# def sample_by_size(adata, props, seed, output_column):
#     # Set seed for reproducibility
#     np.random.seed(seed)  
    
#     # List to store the sampled subsets
#     samples = {}
    
#     # Ensure each class has at least two samples
#     class_groups = adata.obs.groupby(output_column, observed=False)
    
#     for prop in props:
#         # Name for the subset
#         sample_name = f"proportion_{str(prop).replace('.', '_')}"
        
#         # Calculate the total number of samples for this proportion
#         total_size = int(len(adata) * prop)
        
#         # Initialize storage for indices
#         selected_indices = []
        
#         # First, ensure at least 2 samples per class
#         for class_name, class_data in class_groups:
#             if len(class_data) >= 2:
#                 selected_indices += class_data.sample(n=2, random_state=seed).index.tolist()
#             else:
#                 raise ValueError(f"Class {class_name} has fewer than 2 samples, cannot guarantee at least 2 per class.")
        
#         # Remaining size to sample after ensuring 2 samples per class
#         remaining_size = total_size - len(selected_indices)
        
#         if remaining_size > 0:
#             # Create a DataFrame excluding the already selected indices
#             remaining_data = adata[~adata.obs.index.isin(selected_indices)]
            
#             # Randomly sample remaining indices
#             additional_indices = remaining_data.obs.sample(n=remaining_size, random_state=seed).index.tolist()
#             selected_indices += additional_indices
        
#         # Store the sampled subset
#         samples[sample_name] = adata[selected_indices]
    
#     return samples




# -----------------------------------------------------------------------------------------------
# Unused functions
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

# Calculate available CPUs on a host
def get_available_cpus(hostname):
    """
    Get the number of available CPUs for a given host.

    Parameters:
        hostname (str): The name of the host to query.
        save_cpus (int): Number of CPUs to reserve (default is 0).

    Returns:
        int: The number of available CPUs.
    """
    try:
        # Get the qhost output
        qhost_output = subprocess.run(
            ["qhost", "-h", hostname],
            check=True,
            text=True,
            capture_output=True
        ).stdout

        # Match relevant information from the output
        pattern = re.compile(r'(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+(\S+)')
        matches = pattern.findall(qhost_output)

        if not matches:
            raise ValueError(f"No matching host information found for {hostname}")

        # Extract values
        host, arch, ncpu, nsoc, ncor, nthr, load, rest = matches[0]

        # Calculate available CPUs
        ncpu = float(ncpu)
        available_cpus = ncpu

        return available_cpus

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running qhost: {e.stderr.strip()}") from e
    except Exception as e:
        raise RuntimeError(f"An error occurred: {str(e)}") from e


