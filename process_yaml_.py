import yaml
import sys
import os

def modify_yaml(input_path, output_path, seed, njobs):
    """
    Modifies the settings.yml file to submit more jobs to the queuing system.
    """
    with open(input_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Modify the seed and njobs
    data['seed'] = seed
    data['njobs'] = njobs

    # Modify output directory
    base_path = 'data/runs'
    if data['tag'] == 'dev':
        modified_path = 'Dev'

        # Add path to the settings
        data['output_directory'] = os.path.join(base_path, modified_path)
    else:
        # Create a new direcory for each seed
        tag = data['tag']
        comparison = "_vs_".join(data['classification']['comparison'])
        train_ancestry = data['classification']['train_ancestry'].upper()
        infer_ancestry = data['classification']['infer_ancestry'].upper()
        modified_path = f"{tag}_{comparison}_{train_ancestry}_to_{infer_ancestry}_{seed}"

        # Add path to the settings
        data['output_directory'] = os.path.join(base_path, modified_path)

    # Save temp settings file
    with open(output_path, 'w') as f:
        yaml.dump(data, f)

    # Print the modified path to use as qsub output
    print(modified_path)

# Here starts the script

# Seed and input yaml from command line
input_yaml_path = sys.argv[1]  
seed = int(sys.argv[2])
njobs = int(sys.argv[3])         
output_yaml_path = f'_job_settings_{seed}_.yml'  

modify_yaml(input_yaml_path, output_yaml_path, seed, njobs)