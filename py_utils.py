# Standard libraries
import numpy as np
import pandas as pd
import yaml
import re
import math

import uuid
import time
from pathlib import Path
import psutil
import os
import subprocess
import pickle
import importlib
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score

# Data library
import anndata as ad

# Object that contains settings and methods
class Setup(dict):
    def __init__(self, config_file):
        
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
            assert i in final_config, f'No {i} defined but required!'     
        
        # Check required settings in classification
        required_settings = ['comparison', 'output_column', 'train_ancestry', 'infer_ancestry', 'ancestry_column']
        for i in required_settings:
             assert i in final_config['classification'], f'No {i} in classification but required!'
        
        # Check that covariate has to be a list
        assert "covariate" not in final_config["classification"] or isinstance(final_config["classification"]["covariate"], list), \
            "If covariate exists, it must be a list."

        # Check for classification task
        if len(final_config['classification']['comparison']) > 2:
            final_config['classification'].update({'multiclass': True})
            print(f'Multiclass problem comparing: {final_config['classification']['comparison']}.')

        elif len(final_config['classification']['comparison']) == 2:
            final_config['classification'].update({'multiclass': False})
            print(f'Binaryclass problem comparing: {final_config['classification']['comparison']}.')
        
        else:
            print(f'Cant do a classification with one assigned class in comparison: {final_config['classification']['comparison']}.')

        # Init
        self.final_config = final_config
        super().__init__(final_config)

        # Check all the settings
        self.check_settings()

        # Create output directory
        Path(final_config['output_directory']).mkdir(parents=True, exist_ok=self.overwrite)
        save_location = os.path.join(os.getcwd(), final_config['output_directory'])
        # Print tatement to inform about save location
        print(f'Output will be saved to {save_location}')

        # Create a log file
        with open(self.out("Log.tsv"),"w") as file:
            file.write("Step\tMemory_MB\tTime\n")
        
    def __getattr__(self, name):
        return self[name]
    
    def out(self, x):
        return os.path.join(self.output_directory, x)
    
    def check_settings(self):
        # TODO - Check whether settings are of the right type 
        assert isinstance(self.seed, int)

        # Check data_path
        assert isinstance(self.data_path, str)
        assert os.path.exists(self.data_path)
        assert self.data_path.endswith('.h5ad'), f'Only support data files in h5ad format.'

        # Classification
        assert isinstance(self.classification['train_ancestry'], str)
        assert isinstance(self.classification['infer_ancestry'], str)
        assert isinstance(self.classification['ancestry_column'], str)

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

class DataValidator():
    """
    Checks if data and setttings are compatible.
        
    :param data: AnnData object containing the data.
    :param setup: Setup object containing settings.
    """
    def __init__(self, data: ad.AnnData, setup: Setup):

        self.data = data
        self.setup = setup

    def validate_features(self):
        """
        Check that features are unique.
        """
        # Ensure that the number of unique features is equal to the total number of features
        assert self.data.var.index.shape[0] == len(set(self.data.var.index)), \
            'Features not unique.'
    
    def validate_output_column(self):
        """
        Check that output column is in the data.
        """
        output_column = self.setup['classification']['output_column']
        assert output_column in self.data.obs.columns, \
            f"Output column '{output_column}' not in the data."
        
        # Check if there are any NA values 
        has_na = self.data.obs[output_column].isnull().any()
        assert not has_na, \
            f"The covariate column '{output_column}' contains NA values. Please handle missing data before proceeding."

    def validate_class_labels(self):
        """
        Check that class labels of comparisons are present in the output column of the data.
        """
        required_labels = self.setup['classification']['comparison']
        output_column = self.setup['classification']['output_column']
        
        # Ensure that each comparison label is found in the output column
        for label in required_labels:
            assert self.data.obs[output_column].str.contains(label).any(), \
                f"No '{label}' in output column: '{output_column}', choose different comparison in settings."
            
    def validate_ancestry_column(self):
        """
        Check if the ancestry column is present in the data.
        """
        ancestry_column = self.setup['classification']['ancestry_column']
        assert ancestry_column in self.data.obs.columns, \
            f"Ancestry column '{ancestry_column}' not in the data."
        
        # Check if there are any NA values 
        has_na = self.data.obs[ancestry_column].isnull().any()
        assert not has_na, \
            f"The covariate column '{ancestry_column}' contains NA values. Please handle missing data before proceeding."
    
    def validate_ancestry(self):
        """
        Check if ancestries are present in ancetsry column.
        """
        ancestry_column = self.setup['classification']['ancestry_column']
        train_ancestry = self.setup['classification']['train_ancestry']
        inf_ancestry = self.setup['classification']['infer_ancestry']

        assert self.data.obs[ancestry_column].str.contains(train_ancestry).any(), \
            f"No '{train_ancestry}' in ancestry_column: '{ancestry_column}'."
        assert self.data.obs[ancestry_column].str.contains(inf_ancestry).any(), \
            f"No '{inf_ancestry}' in ancestry_column: '{ancestry_column}'."

    def validate_covariate_column(self):
        """
        Check if covariate column(s) are present in the data, if specified.
        Check that there are no NA values in the covariate column(s).
        """
        classification_settings = self.setup.get('classification', {})
        
        # Check if 'covariate' is provided in the classification settings
        covariates = classification_settings.get('covariate')
        if covariates:  # If covariates are specified
            if not isinstance(covariates, list):
                covariates = [covariates]  # Ensure it's a list
            
            for covariate in covariates:
                # Check if each covariate column exists in the data
                assert covariate in self.data.obs.columns, \
                    f"No '{covariate}' column in the data, can't be used as covariate."
                
                # Check if there are any NA values in the covariate column
                has_na = self.data.obs[covariate].isnull().any()
                assert not has_na, \
                    f"The covariate column '{covariate}' contains NA values. Please handle missing data before proceeding."

        
    # Function that runs all validate functions
    def validate(self):
        """
        Run all validation checks on the data.
        """
        self.validate_features()
        self.validate_output_column()
        self.validate_class_labels()
        self.validate_ancestry_column()
        self.validate_ancestry()
        self.validate_covariate_column()


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

        # Save the best hyperparameters to a CSV file
        hyp = pd.DataFrame(grid_search.best_params_, index=[setup["id"]])
        if prop is not None:
            hyp.to_csv(os.path.join(setup["output_directory"], f"Hyperparameters_{prop}_{algo_name}.csv"), index=False)
            # Optionally save the trained model to a pickle file
            if setup["save_model"]:
                with open(os.path.join(setup["output_directory"], f"{algo_name}_{prop}_{setup["id"]}.pkl"), "wb") as f:
                    pickle.dump(best_model, f, protocol=pickle.HIGHEST_PROTOCOL)

        else:
            hyp.to_csv(os.path.join(setup["output_directory"], f"Hyperparameters_{algo_name}.csv"), index=False)
            # Optionally save the trained model to a pickle file
            if setup["save_model"]:
                with open(os.path.join(setup["output_directory"], f"{algo_name}_{setup["id"]}.pkl"), "wb") as f:
                    pickle.dump(best_model, f, protocol=pickle.HIGHEST_PROTOCOL)

        # Training complete
        print(f"{algo_name} training done.")
        return best_model

    except Exception as e:
        # Handle and log any exceptions that occur during training
        print(f"Error occurred during training for {algo_name}: {e}")
        return None
    

def evaluate_algorithm_cross_ancestry(algo_name, best_model, setup, test_X, test_y, inf_X, inf_y, feature_names):
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

    Returns:
        None
    """
    try:
        # Log the start of testing
        setup.log("Testing")

        # Generate predictions on the test set and save probabilities
        test_y_hat = pd.DataFrame(best_model.predict_proba(test_X))
        test_y_hat["y"] = test_y
        test_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_test.csv"), index=False)

        # Log the start of inference
        setup.log("Infering")

        # Generate predictions on the inference set and save probabilities
        inf_y_hat = pd.DataFrame(best_model.predict_proba(inf_X))
        inf_y_hat["y"] = inf_y
        inf_y_hat.to_csv(setup.out(f"Probabilities_{algo_name}_inf.csv"), index=False)

        # Inference complete
        print(f"{algo_name} validation done.")

        # Log the start of model interpretation
        setup.log("Model interpretations")

        # Extract feature importances (if supported by the model)
        feature_importance = extract_feature_importance(best_model, feature_names)

        # Save feature importances to a CSV file
        feature_importance.to_csv(setup.out(f"Feature_importance_{algo_name}.csv"), index=False)

    except Exception as e:
        # Handle and log any exceptions that occur during evaluation
        print(f"Error occurred during evaluation for {algo_name}: {e}")


def evaluate_algorithm_eur_subsetting(algo_name, setup, best_model, subset_X, subset_y):
    """
    Evaluate the trained model on test subsets and return probabilities and metrics.
    """
    try:
        propability_list = []
        metric_list = []
        
        # Iterate over all subsets (X, y)
        for X, y in zip(subset_X, subset_y):
            # Get predicted probabilities
            test_y_hat = pd.DataFrame(best_model.predict_proba(X))
            test_y_hat["y"] = y
            test_y_hat["n_test_ancestry"] = y.shape[0]
            
            # Collect probabilities
            propability_list.append(test_y_hat)

            # Calculate metrics
            test_probabilities, test_y = test_y_hat.iloc[:, 1].values, test_y_hat["y"].values
            auc = roc_auc_score(y_true=test_y, y_score=test_probabilities)
            metric = {algo_name: auc, "n_test_ancestry": y.shape[0], "Metric": "ROC_AUC"}
            
            # Collect metrics
            metric_list.append(metric)

        # Combine probabilities and save
        propabilities = pd.concat(propability_list, ignore_index=True)
        propabilities.to_csv(setup.out(f"Probabilities_{algo_name}.csv"), index=False)

        # Combine metric and add algo information
        metric_combined = pd.DataFrame(metric_list)

        # Finished validation
        print(f"{algo_name} validation done.")

        return metric_combined
    
    except Exception as e:
        print(f"Error occurred during evaluation for {algo_name}: {e}")



def evaluate_algorithm_robustness(algo_name, best_model, setup, test_X, test_y, inf_X, inf_y, feature_names, prop):
    """
    Evaluate the trained model, save probabilities, calculate metrics, and return results.
    """
    try:
        setup.log("Testing")
        # Testing on the test set
        test_y_hat = pd.DataFrame(best_model.predict_proba(test_X))
        test_y_hat["y"] = test_y
        test_y_hat.to_csv(setup.out(f"Probabilities_{prop}_{algo_name}_test.csv"), index=False)

        setup.log("Infering")
        # Testing on the inference set
        inf_y_hat = pd.DataFrame(best_model.predict_proba(inf_X))
        inf_y_hat["y"] = inf_y
        inf_y_hat.to_csv(setup.out(f"Probabilities_{prop}_{algo_name}_inf.csv"), index=False)

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

        return metric_df

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
    importance_df = pd.DataFrame({
        "Feature": feature_names,
        "Importance": importance
    }).sort_values(by="Importance", ascending=False)
    
    return importance_df





# EUR Subsetting / Robustness
def sample_by_size(adata, props, seed, output_column):
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
            else:
                raise ValueError(f"Class {class_name} has fewer than 2 samples, cannot guarantee at least 2 per class.")
        
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


