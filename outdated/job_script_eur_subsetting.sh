#!/bin/bash

# Here starts the script

# Make the settings and seeds command line argument
# First argument
SETTINGS_FILE="$1"
shift
# List of seeds to run (== amount of jobs)
# Second argument
seeds=("$@")

# Need to specify host
HOSTNAME="trude.came.sbg.ac.at"
SAVE_CPUS=0

# Process settings file
PROCESS_SCRIPT="process_yaml_eur_subsetting_.py"
# Singularity container for execution
SINGULARITY_IMAGE="data/ancestry_dk.sif"

# Scripts to run
SCRIPT_SUBSETTING="main_eur_subsetting.py"
R_SCRIPT_COMBINE="combine_runs.R"


# Check if the YAML file exists
if [ ! -f "$SETTINGS_FILE" ]; then
  echo "Error: YAML file does not exist at $SETTINGS_FILE"
  exit 1
fi

# Make scripts executable
chmod +x "${PROCESS_SCRIPT}"

chmod +x "${SCRIPT_SUBSETTING}"
chmod +x "${R_SCRIPT_COMBINE}"


# Check availability of cpus on host
qhost_output=$(qhost -h "${HOSTNAME}")

# Parse the output (skip the header and extract relevant information)
read -r HOST ARCH NCPU NSOC NCOR NTHR LOAD REST <<<$(echo "$qhost_output" | tail -n +4)

# Calculate available CPUs based on load, considering the reserved CPUs
AVAILABLE_CPUS=$(echo "($NCPU - $LOAD) - ${SAVE_CPUS}" | bc)
AVAILABLE_CPUS=$(printf "%.0f" "$AVAILABLE_CPUS")

# Cpus per seed
CPUS_PER_JOB=$((AVAILABLE_CPUS / ${#seeds[@]}))


job_list=()
for seed in "${seeds[@]}"; do
    # echo "$seed"
    # Modify settings file (singularity does not work within container)
    job_name=$(singularity exec --bind /vscratch:/vscratch "${SINGULARITY_IMAGE}" python3 "${PROCESS_SCRIPT}" "${SETTINGS_FILE}" "${seed}" "${CPUS_PER_JOB}")

    job_list+=${job_name},
    # Make output for qsub
    qsub_output="data/qsub/$job_name"

    # Check if the directory exists first (should raise error when already exists)
    if [ ! -d "$qsub_output" ]; then
        mkdir -p "$qsub_output"
    else
        rm -r "$qsub_output"
        mkdir -p "$qsub_output"
    fi

    # Submit the job to the queuing system
    # singularity exec --bind /vscratch:/vscratch "${SINGULARITY_IMAGE}" python3 "${PYTHON_SCRIPT}" job_settings_.yml
    
    # This is the one that runs the actual job
    qsub -N $job_name -l hostname=$HOSTNAME -pe smp $CPUS_PER_JOB -o "$qsub_output/output.log" -e "$qsub_output/error.log" singularity exec --bind /vscratch:/vscratch "${SINGULARITY_IMAGE}" python3 "${SCRIPT_SUBSETTING}" _job_settings_"${seed}_".yml

    # qsub -N $job_name -l hostname=$HOSTNAME -o "$qsub_output/output.log" -e "$qsub_output/error.log" singularity exec --bind /vscratch:/vscratch "${SINGULARITY_IMAGE}" python3 "${PYTHON_SCRIPT}" job_settings_.yml
    # qsub -N "echo" -l hostname=$HOSTNAME -o "$qsub_output/output.log" -e "$qsub_output/error.log" echo "HALLO_there_should_be_no_error"
done

# These scripts wait to be executed once all jobs in job_list are finished
# qsub -N "Combine_runs" -hold_jid "${job_list[@]}" singularity exec --bind /vscratch:/vscratch "${SINGULARITY_IMAGE}" Rscript "${R_SCRIPT}" "${SETTINGS_FILE}"

