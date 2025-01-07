import os
import sys
import stat
import shutil
import subprocess

import yaml
import numpy as np
import re

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


# Settings that should not change 
SINGULARITY = "data/ancestry_dk.sif"
HOSTNAME = "trude.came.sbg.ac.at"

# This script should execute 'robustness analysis'
# Inputs:
#   Settings:   YAML file
#   Seeds:      List of seeds 

SETTINGS_FILE = "EUR_ASI_settings.yml"
SEEDS = np.arange(1,11)

# Check if 'SETTINGS_FILE' exists
if not os.path.isfile(SETTINGS_FILE):
    print(f"Error: YAML file does not exist at {SETTINGS_FILE}")
    sys.exit(1)
# Load the 'SETTINGS_FILE'
else:
    with open(SETTINGS_FILE, 'r') as f:
        settings = yaml.safe_load(f)

# Dynamically assess CPUs on 'HOSTNAME'
# Calculate CPUs on hostname (save some cpus)
save_cpus = 10
available_cpus = get_available_cpus(hostname=HOSTNAME) - save_cpus


# Calculate CPUs per job
cpus_per_job = int(max(1, available_cpus // len(SEEDS)))
settings['n_jobs'] = cpus_per_job

# Run 'robustness' analysis
script_name = "robustness.py"

# Define analaysis name
train_ancestry = settings['classification']['train_ancestry'].upper()
infer_ancestry = settings['classification']['infer_ancestry'].upper()
analysis_name = f"{train_ancestry}_to_{infer_ancestry}_robustness"

# Make script executable
os.chmod(script_name, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

# Save job names
job_name_list = list()
for seed in SEEDS:
    # Modify seed
    settings["seed"] = int(seed)

    # Modify 'output_directory'
    # 'data/runs' is the directory where individual seeds are saved
    vscratch_dir = os.path.join("data", "runs")
    # Create directory for each seed
    # Tag is used to specify which data it is e.g. TCGA, NIAGADS
    tag = settings['tag']
    # Comparison is specified which conditions are compared e.g. cancer_1_vs_cancer_2
    comparison = "_vs_".join(settings['classification']['comparison'])

    # 'analysis_name' per seed
    analysis_name_seed = f"{analysis_name}_{seed}"

    # Final modified path
    # 'modified_path' is used as job_name
    modified_path = f"{tag}_{comparison}_{analysis_name_seed}"
    # Add modified path to settings
    settings['output_directory'] = os.path.join(vscratch_dir, modified_path)
    # Add 'modified_path' to 'job_name_list'
    job_name_list.append(modified_path)

    # Create 'tmp' directory on vsctrach
    directory_name = os.path.join("data", "tmp")
    os.makedirs(directory_name, exist_ok=True)

    # Save the settings file so it can be used for qsub
    path_to_modified_settings = os.path.join(directory_name, f"_{modified_path}_.yml")
    with open(path_to_modified_settings, "w") as f:
        yaml.dump(settings, f)
    
    # Create 'qsub' directory on vscratch
    qsub_directory = os.path.join("data", "qsub")
    os.makedirs(qsub_directory, exist_ok=True)

    # Save qsub output
    qsub_output = os.path.join(qsub_directory, modified_path)
    if not os.path.isdir(qsub_output):
        os.makedirs(qsub_output)
    else:
        # If it exists, remove it and recreate it
        shutil.rmtree(qsub_output)  
        os.makedirs(qsub_output)  

    # Define the qsub command
    command = [
    "qsub",
    "-N", modified_path,
    "-l", f"hostname={HOSTNAME}",
    "-pe", "smp", str(cpus_per_job),
    "-o", f"{qsub_output}/output.log",
    "-e", f"{qsub_output}/error.log",
    "singularity", "exec", "--bind", "/vscratch:/vscratch",
    SINGULARITY,
    "python3", script_name,
    path_to_modified_settings]

    # Execute the command
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"Command executed successfully:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command:\n{e.stderr}")
