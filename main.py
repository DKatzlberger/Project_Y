import os
import sys
import stat
import shutil
import subprocess
import yaml
import numpy as np
import re
import time

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

def submit_seed_jobs(settings_file, 
                     seeds, 
                     hostname, 
                     singularity, 
                     script_name, 
                     analysis_name, 
                     cpus=10):
    """
    Prepare and submit seed jobs for analysis.

    Parameters:
        settings_file (str): Path to the settings YAML file.
        seeds (list): List of seeds for the analysis.
        hostname (str): Hostname for CPU availability.
        singularity (str): Path to the Singularity image.
        script_name (str): Name of the script to execute.
        save_cpus (int): Number of CPUs to reserve.
        analysis_name (str): Optional dynamic name for the analysis. If None, it will be generated dynamically.

    Returns:
        list: A list of job names corresponding to the submitted seed jobs.
    """
    # Load settings
    settings = load_settings(settings_file)

    # Set cpus
    settings["njobs"] = cpus

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
            "-pe", "smp", str(cpus),
            "-o", os.path.join(qsub_output_dir, "output.log"),
            "-e", os.path.join(qsub_output_dir, "error.log"),
            "singularity", "exec", "--bind", "/vscratch:/vscratch",
            singularity,
            "python3", script_name, modified_settings_path
        ]

        # Execute the command
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            print(f"Command executed successfully:\n{result.stdout}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command:\n{e.stderr}")

    return job_name_list, modified_settings_path


def submit_single_job(settings_file,
                      script_name, 
                      analysis_name, 
                      singularity,
                      cpus=1,
                      hostname=None,
                      job_name_list=None,
                      wait_for_jobs=True
                      ):
    """
    Submit a job that optionally waits for seed jobs to finish before combining their outputs.

    Parameters:
        settings_file (str): Path to the settings file.
        script_name (str): Name of the script to execute.
        analysis_name (str): Name of the analysis.
        singularity (str): Path to the Singularity image.
        cpus (int): Number of CPU cores requested.
        hostname (str, optional): Hostname to run the job on.
        job_name_list (list, optional): List of job names to wait for.
        wait_for_jobs (bool): Whether to wait for the jobs in job_name_list to finish.

    Returns:
        str: The name of the submitted job.
    """

    # Load settings
    settings = load_settings(settings_file)

    # Create job name
    tag = settings['tag']
    comparison = "_vs_".join(settings['classification']['comparison'])
    job_name = f"{tag}_{comparison}_{analysis_name}"

    # Define paths and directories
    tmp_dir = os.path.join("data", "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    qsub_dir = os.path.join("data", "qsub")
    os.makedirs(qsub_dir, exist_ok=True)

    # Prepare qsub output directory
    qsub_output_dir = prepare_directories(qsub_dir, job_name)

    # Make script executable
    os.chmod(script_name, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

    # Construct qsub command
    command = [
        "qsub",
        "-N", job_name,
        "-o", os.path.join(qsub_output_dir, "output.log"),
        "-e", os.path.join(qsub_output_dir, "error.log"),
        "-pe", "smp", str(cpus),  # Request CPUs
        "singularity", "exec", "--bind", "/vscratch:/vscratch",
        singularity,
        "Rscript", script_name, settings_file
    ]

    # Add hostname constraint if provided
    if hostname:
        command.insert(3, "-l")
        command.insert(4, f"hostname={hostname}")

    # Add hold_jid only if wait_for_jobs is True and job_name_list is not empty
    if wait_for_jobs and job_name_list:
        ids_to_wait = ",".join(job_name_list)
        command.insert(3, "-hold_jid")
        command.insert(4, ids_to_wait)

    # Execute the command
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"Command executed successfully:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command:\n{e.stderr}")

    return job_name



# Main execution
if __name__ == "__main__":

    # Fixed arguments
    SINGULARITY = "data/ancestry_dk.sif"
    HOSTNAME = "trude.came.sbg.ac.at"
    SEEDS = np.arange(1, 11)

    # Ancestry specific settings
    SETTINGS_FILES = [
        "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_ADMIX.yml",
        "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_AFR.yml", 
        "data/inputs/settings/PanCanAtlas_BRCA_RSEM_basal_vs_non-basal_EUR_to_EAS.yml"
        # "data/inputs/settings/PanCanAtlas_LUSC_LUAD_RSEM_Lung_Adenocarcinoma_vs_Lung_Squamous_Cell_Carcinoma_EUR_to_ADMIX.yml",
        # "data/inputs/settings/PanCanAtlas_LUSC_LUAD_RSEM_Lung_Adenocarcinoma_vs_Lung_Squamous_Cell_Carcinoma_EUR_to_AFR.yml",
        # "data/inputs/settings/PanCanAtlas_LUSC_LUAD_RSEM_Lung_Adenocarcinoma_vs_Lung_Squamous_Cell_Carcinoma_EUR_to_EAS.yml"
    ]
    # Pattern to extract 'eur_to_region'
    pattern = r"EUR_to_[A-Z]+"

    # CPUS
    available_cpus = get_available_cpus(hostname=HOSTNAME, save_cpus=5)
    number_jobs = len(SETTINGS_FILES) * len(SEEDS) * 3
    available_cpus = max(4, available_cpus // number_jobs)

    # Job list for all ancestries
    # 'Combine_ancestries' script needs to wait for these jobs
    combine_ancestries_wait_jobs = []
    for SETTINGS_FILE in SETTINGS_FILES:
        
        # Extract ancestry comparison
        match = re.search(pattern, SETTINGS_FILE)
        eur_to_region = match.group()

        # 'cross_ancestry' analysis
        SCRIPT = "cross_ancestry.py"
        ANALYSIS_NAME = f"{eur_to_region}_cross_ancestry"

        # Submit seed jobs and retrieve the job names
        cross_ancestry_jobs, tmp_settings = submit_seed_jobs(
            settings_file=SETTINGS_FILE,
            seeds=SEEDS,
            hostname=HOSTNAME,
            singularity=SINGULARITY,
            script_name=SCRIPT,
            analysis_name=ANALYSIS_NAME,
            cpus=available_cpus
        )

        # 'Interactions'
        INTERACTIONS_SCRIPT = "interactions.R"
        NAME = f"{eur_to_region}_interactions"

        interactions = submit_single_job(
            settings_file=tmp_settings,
            script_name=INTERACTIONS_SCRIPT,
            analysis_name=NAME,
            singularity=SINGULARITY,
            hostname=HOSTNAME,
            wait_for_jobs=True,
            job_name_list=cross_ancestry_jobs,
            cpus=available_cpus
        )

        # Combine wait jobs
        interactions_wait_jobs = cross_ancestry_jobs
        interactions_wait_jobs.append(interactions)

        # 'cross_ancestry' combine runs
        COMBINE_SCRIPT = "cross_ancestry_combine_runs_01.R"
        NAME = f"{eur_to_region}_cross_ancestry_combine_runs"

        cross_ancestry_combine_runs_job = submit_single_job(
            settings_file=tmp_settings,
            script_name=COMBINE_SCRIPT,
            analysis_name=NAME,
            singularity=SINGULARITY,
            hostname=HOSTNAME,
            wait_for_jobs=True,
            job_name_list=interactions_wait_jobs,
            cpus=available_cpus
        )

        # 'robustness' analysis
        SCRIPT = "robustness.py"
        ANALYSIS_NAME = f"{eur_to_region}_robustness"

        # Submit seed jobs and retrieve the job names
        robustness_jobs, tmp_settings = submit_seed_jobs(
            settings_file=SETTINGS_FILE,
            seeds=SEEDS,
            hostname=HOSTNAME,
            singularity=SINGULARITY,
            script_name=SCRIPT,
            analysis_name=ANALYSIS_NAME,
            cpus=available_cpus
        )

        # 'robustness' combine runs
        COMBINE_SCRIPT = "robustness_combine_runs_01.R"
        NAME = f"{eur_to_region}_robustness_combine_runs"

        robustness_combine_runs_job = submit_single_job(
            settings_file=tmp_settings,
            script_name=COMBINE_SCRIPT,
            analysis_name=NAME,
            singularity=SINGULARITY,
            hostname=HOSTNAME,
            wait_for_jobs=True,
            job_name_list=robustness_jobs,
            cpus=1
        )

        # Add to combine_ancestries_wait_jobs
        combine_ancestries_wait_jobs.extend([cross_ancestry_combine_runs_job, robustness_combine_runs_job])


    # This analysis only needs to run once per comparison
    # "descriptive_model_building"
    SCRIPT_NAME = "descriptive_model_building.R"
    NAME = "descriptive_statisitics"
    job_name = submit_single_job(
            settings_file=SETTINGS_FILE,
            script_name=SCRIPT_NAME,
            analysis_name=NAME,
            singularity=SINGULARITY,
            hostname=HOSTNAME,
            wait_for_jobs=False,
            cpus=1
        )


    # 'eur_subsetting' analysis
    SCRIPT_NAME = "eur_subsetting.py"
    NAME = "EUR_subsetting"
    # Submit seed jobs and retrieve the job names
    eur_subsetting_jobs, tmp_settings = submit_seed_jobs(
        settings_file=SETTINGS_FILE,
        seeds=SEEDS,
        hostname=HOSTNAME,
        singularity=SINGULARITY,
        script_name=SCRIPT_NAME,
        analysis_name=NAME,
        cpus=available_cpus
    )

    # 'eur-subsetting' combine runs
    COMBINE_SCRIPT_NAME = "eur_subsetting_combine_runs.R"
    NAME = f"{eur_to_region}_eur_subsetting_combine_runs"

    # Submit combine job that depends on seed jobs
    eur_subsetting_combine_runs_job = submit_single_job(
        settings_file=tmp_settings,
        script_name=COMBINE_SCRIPT_NAME,
        analysis_name=NAME,
        singularity=SINGULARITY,
        hostname=HOSTNAME,
        wait_for_jobs=True,
        job_name_list=eur_subsetting_jobs,
        cpus=1
    )

    # 'Combine ancestries' script
    COMBINE_ANCESTRIES_SCRIPT = "cross_ancestry_combine_ancestries.R"
    NAME = "combine_ancestries"
    combine_ancestries_wait_jobs.append(eur_subsetting_combine_runs_job)
    
    job_name = submit_single_job(
            settings_file=SETTINGS_FILE,
            script_name=COMBINE_ANCESTRIES_SCRIPT,
            analysis_name=NAME,
            singularity=SINGULARITY,
            hostname=HOSTNAME,
            wait_for_jobs=True,
            job_name_list=combine_ancestries_wait_jobs,
            cpus=1
        )
    