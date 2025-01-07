import os
import sys
import stat
import shutil
import subprocess
import yaml
import numpy as np
import re

def get_available_cpus(hostname, save_cpus=10):
    """
    Get the number of available CPUs for a given host, minus reserved CPUs.

    Parameters:
        hostname (str): The name of the host to query.
        save_cpus (int): Number of CPUs to reserve.

    Returns:
        int: The number of available CPUs.
    """
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

    _, _, ncpu, _, _, _, _, _ = matches[0]

    # Calculate available CPUs
    ncpu = int(ncpu)
    available_cpus = max(0, int(ncpu - save_cpus))
    return int(available_cpus)


def load_settings(file_path):
    """
    Load a YAML settings file.

    Parameters:
        file_path (str): Path to the YAML file.

    Returns:
        dict: The settings loaded from the YAML file.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Error: YAML file does not exist at {file_path}")
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)
    

def save_settings(settings, file_path):
    """
    Save a settings dictionary to a YAML file.

    Parameters:
        settings (dict): The settings to save.
        file_path (str): Path to the YAML file.
    """
    with open(file_path, "w") as f:
        yaml.dump(settings, f)


def make_executable(file_path):
    """
    Make a script executable by granting all permissions.

    Parameters:
        file_path (str): Path to the file to make executable.
    """
    os.chmod(file_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)


def prepare_directories(base_dir, sub_dir):
    """
    Create a directory, removing it first if it exists.

    Parameters:
        base_dir (str): Base directory path.
        sub_dir (str): Subdirectory to create.
    """
    directory = os.path.join(base_dir, sub_dir)
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)
    return directory

def prepare_and_submit_jobs(settings_file, seeds, hostname, singularity, script_name, save_cpus=10, analysis_name=None):
    """
    Prepare and submit jobs for robustness analysis with dynamic analysis name.

    Parameters:
        settings_file (str): Path to the settings YAML file.
        seeds (list): List of seeds for the analysis.
        hostname (str): Hostname for CPU availability.
        singularity (str): Path to the Singularity image.
        script_name (str): Name of the script to execute.
        save_cpus (int): Number of CPUs to reserve.
        analysis_name (str): Optional dynamic name for the analysis. If None, it will be generated dynamically.
    """
    # Load settings
    settings = load_settings(settings_file)

    # Calculate available CPUs
    available_cpus = get_available_cpus(hostname, save_cpus)
    cpus_per_job = max(1, available_cpus // len(seeds))
    settings['n_jobs'] = cpus_per_job

    # Define paths and constants
    tmp_dir = os.path.join("data", "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    qsub_dir = os.path.join("data", "qsub")
    os.makedirs(qsub_dir, exist_ok=True)
    
    # Make script executable
    make_executable(script_name)

    job_name_list = []
    for seed in seeds:
        # Update settings
        settings["seed"] = int(seed)
        tag = settings['tag']
        comparison = "_vs_".join(settings['classification']['comparison'])
        analysis_name_seed = f"{analysis_name}_{seed}"
        modified_path = f"{tag}_{comparison}_{analysis_name_seed}"
        settings['output_directory'] = os.path.join("data", "runs", modified_path)
        job_name_list.append(modified_path)

        # Save modified settings
        modified_settings_path = os.path.join(tmp_dir, f"_{modified_path}_.yml")
        save_settings(settings, modified_settings_path)

        # Prepare qsub output directory
        qsub_output_dir = prepare_directories(qsub_dir, modified_path)

        # Submit job using qsub
        command = [
            "qsub",
            "-N", modified_path,
            "-l", f"hostname={hostname}",
            "-pe", "smp", str(cpus_per_job),
            "-o", os.path.join(qsub_output_dir, "output.log"),
            "-e", os.path.join(qsub_output_dir, "error.log"),
            "singularity", "exec", "--bind", "/vscratch:/vscratch",
            singularity,
            "python3", script_name,
            modified_settings_path
        ]

        # Execute the command
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            print(f"Command executed successfully:\n{result.stdout}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command:\n{e.stderr}")



# Main execution
if __name__ == "__main__":
    SINGULARITY = "data/ancestry_dk.sif"
    HOSTNAME = "trude.came.sbg.ac.at"
    SEEDS = np.arange(1, 11)
    SETTINGS_FILE = "EUR_to_ASI_settings.yml"
    SCRIPT_NAME = "robustness.py"

    # You can define a custom analysis name or leave it to be generated dynamically
    CUSTOM_ANALYSIS_NAME = "EUR_to_ASI_robustness"

    prepare_and_submit_jobs(
        SETTINGS_FILE,
        SEEDS,
        HOSTNAME,
        SINGULARITY,
        SCRIPT_NAME,
        save_cpus=10,
        analysis_name=CUSTOM_ANALYSIS_NAME
    )