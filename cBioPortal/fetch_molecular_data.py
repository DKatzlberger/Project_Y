# This script needs a working conda environment with pybioportal installed
# Pybioportal
# import pybioportal 
from pybioportal import molecular_data as md
from pybioportal import molecular_profiles as mpf
from pybioportal import server_running_status as srs
# Standard libraries
import pandas as pd
import os

# Fetches molecular data from the cBioPortal API
# Meta data that contains study_ids and sample_ids
# 'cBioPortal' files are in a different folder to have access:
# One step up in the folder hierachy
to_fetch = pd.read_csv("../data/downloads/cbioportal/tcga_firehose/meta_tcga_firehose_methylation.csv")

# to_fetch = pd.read_csv('data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv')
# study_id = to_fetch['studyid'].unique()[2]

molecular_data = []
for study_id in to_fetch['studyId'].unique():
    #study_id = "brca_tcga_pan_can_atlas_2018"
    print(f"Fetching {study_id}")
    # Molecular profile ids
    # 1. Get all molecular profiles in study
    # 2. Filter for 'mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)'
    profiles = mpf.get_all_molecular_profiles_in_study(study_id=study_id)

    # Available profiles
# 0                      uvm_tcga_pan_can_atlas_2018_rppa
# 1              uvm_tcga_pan_can_atlas_2018_rppa_Zscores
# 2                    uvm_tcga_pan_can_atlas_2018_gistic
# 3                   uvm_tcga_pan_can_atlas_2018_log2CNA
# 4              uvm_tcga_pan_can_atlas_2018_armlevel_cna
# 5                 uvm_tcga_pan_can_atlas_2018_mutations
# 6       uvm_tcga_pan_can_atlas_2018_structural_variants
# 7     uvm_tcga_pan_can_atlas_2018_methylation_hm27_h...
# 8           uvm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna
# 9     uvm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_me...
# 10    uvm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_me...
# 11         uvm_tcga_pan_can_atlas_2018_genetic_ancestry

    met_id = f"{study_id}_methylation_hm450"
    rna_id = f'{study_id}_rna_seq_v2_mrna'
    rna_profile = profiles[profiles['molecularProfileId'] == met_id]
    # Get 'molecularProfileId'
    rna_profile_id = rna_profile['molecularProfileId'].values[0]

    # Samples in the study
    # Used to fetch by individual sample
    to_fetch_samples = to_fetch[to_fetch['studyId'] == study_id]['sampleId'].tolist()
    print(len(to_fetch_samples))

    # Save location
    vscratch_dir_out = f"../data/downloads/cbioportal/tcga_firehose/methylation"
    # 1. Location for successful samples
    # 2. Location for failed samples
    path_to_successful = os.path.join(vscratch_dir_out, f"{study_id}.csv")
    path_to_failed = os.path.join(vscratch_dir_out, f"failed_samples_{study_id}.txt")

    if len(to_fetch_samples) < 500:
        try: 
            # Fetch all molecular data for 'molecularProfileId' and available RNA data
            print(f"Trying to fetch by 'sample_list_id'.")
            response = md.fetch_all_molecular_data_in_molecular_profile(molecular_profile_id=rna_profile_id, 
                                                                        sample_list_id=rna_profile_id,
                                                                        projection='DETAILED')

            # Filter for samples that are in the meta data
            # Because some dont have genetic ancestry
            # 1. Get samples with ancestry from meta data
            # 2. Filter for samples with ancestry
            # 3. Filter for important attributes
            samples_with_ancestry = to_fetch[to_fetch['studyId'] == study_id]['sampleId']
            attributes = ['gene_entrezGeneId', 'gene_hugoGeneSymbol', 'sampleId', 'patientId', 'value', 'studyId']
            # Filter
            response = response[response['sampleId'].isin(samples_with_ancestry)]
            response = response[attributes]

            # Pivot dataframe wider (will drop Entrez gene id)
            response_wide = response.pivot(index=['sampleId', 'patientId', 'studyId'], 
                                        columns='gene_hugoGeneSymbol', 
                                        values='value').reset_index()
            
            # Check if all samples have been downloaded
            if not response_wide['sampleId'].isin(to_fetch_samples).all():
                print(f"Not all samples downloaded.")
            else:
                print("Successful.\n")

            # Save individual study
            response_wide.to_csv(path_to_successful, index=False)

            # Create a list of dataframes
            molecular_data.append(response_wide)

        except:
            # Fetch data by individual samples
            print(f"Trying to fetch by individual sample.")

            # Initialize data frame to combine all samples
            response = pd.DataFrame()
            failed_samples = []
            sample_count = 1
            # Iterate over each sample
            for i in range(0, len(to_fetch_samples)): 

                # 1. Fetch each individual sample
                # 2. Combine them to complete data frame
                # 3. Failed samples append to list
                sample = to_fetch_samples[i:i + 1]
                try:
                    sample_response = md.fetch_all_molecular_data_in_molecular_profile(molecular_profile_id=rna_profile_id,
                                                                                    sample_ids=sample,
                                                                                    projection='DETAILED')
                    response = pd.concat([response, sample_response], axis=0)
                    # Statement to keep track of samples 
                    print(f"{sample_count}/{len(to_fetch_samples)}")

                except:
                    # Append failed sample to 'failed_samples'
                    print(f"'{sample[0]}' failed to download.")
                    failed_samples.append(sample[0])
                
                # Add counter
                sample_count += 1

            # Filter for samples that are in the meta data
            # Because some dont have genetic ancestry
            # 1. Get samples with ancestry from meta data
            # 2. Filter for samples with ancestry
            # 3. Filter for important attributes
            samples_with_ancestry = to_fetch[to_fetch['studyId'] == study_id]['sampleId']
            attributes = ['gene_entrezGeneId', 'gene_hugoGeneSymbol', 'sampleId', 'patientId', 'value', 'studyId']
            # Filter
            response = response[response['sampleId'].isin(samples_with_ancestry)]
            response = response[attributes]

            # Pivot dataframe wider (will drop Entrez gene id)
            response_wide = response.pivot(index=['sampleId', 'patientId', 'studyId'], 
                                        columns='gene_hugoGeneSymbol', 
                                        values='value').reset_index()

            # Check if all samples have been downloaded
            if not response_wide['sampleId'].isin(to_fetch_samples).all():
                print(f"Not all samples downloaded.")
            else:
                print("Successful.\n")

            # Save individual study
            response_wide.to_csv(path_to_successful, index=False)

            # If 'failed_samples' is not empty save them in a list
            if failed_samples:
                with open(path_to_failed, 'w') as file:
                    for item in failed_samples:
                        file.write(item + '\n')

            # Create a list of dataframes
            molecular_data.append(response_wide)

    else:
        # Fetch individual samples if more than 500 samples in study
        # Fetch data by individual samples
        print(f"Study has more than 500 samples.")
        print(f"Trying to fetch by individual sample.")

        # Initialize data frame to combine all samples
        response = pd.DataFrame()
        failed_samples = []
        sample_count = 1
        # Iterate over each sample
        for i in range(0, len(to_fetch_samples)): 

            # 1. Fetch each individual sample
            # 2. Combine them to complete data frame
            # 3. Failed samples append to list
            sample = to_fetch_samples[i:i + 1]
            try:
                sample_response = md.fetch_all_molecular_data_in_molecular_profile(molecular_profile_id=rna_profile_id,
                                                                                sample_ids=sample,
                                                                                projection='DETAILED')
                response = pd.concat([response, sample_response], axis=0)
                # Statement to keep track of samples 
                print(f"{sample_count}/{len(to_fetch_samples)}")

            except:
                # Append failed sample to 'failed_samples'
                print(f"'{sample[0]}' failed to download.")
                failed_samples.append(sample[0])
            
            # Add counter
            sample_count += 1

        # Filter for samples that are in the meta data
        # Because some dont have genetic ancestry
        # 1. Get samples with ancestry from meta data
        # 2. Filter for samples with ancestry
        # 3. Filter for important attributes
        samples_with_ancestry = to_fetch[to_fetch['studyId'] == study_id]['sampleId']
        attributes = ['gene_entrezGeneId', 'gene_hugoGeneSymbol', 'sampleId', 'patientId', 'value', 'studyId']
        # Filter
        response = response[response['sampleId'].isin(samples_with_ancestry)]
        response = response[attributes]

        # Pivot dataframe wider (will drop Entrez gene id)
        response_wide = response.pivot(index=['sampleId', 'patientId', 'studyId'], 
                                    columns='gene_hugoGeneSymbol', 
                                    values='value').reset_index()

        # Check if all samples have been downloaded
        if not response_wide['sampleId'].isin(to_fetch_samples).all():
            print(f"Not all samples downloaded.")
        else:
            print("Successful.\n")

        # Save individual study
        response_wide.to_csv(path_to_successful, index=False)

        # If 'failed_samples' is not empty save them in a list
        if failed_samples:
            with open(path_to_failed, 'w') as file:
                for item in failed_samples:
                    file.write(item + '\n')

        # Create a list of dataframes
        molecular_data.append(response_wide)


print('Finished!')
# Assertion: Check if all studies have been downloaded
if not len(molecular_data) == len(to_fetch['studyId'].unique()):
    print("Not all studies have samples.")

# Combine all dataframes: 
# Using 'pd.concat' would produce NA values in columns (= features) which aren't available across studies
# combined_molecular_data = pd.concat(molecular_data, axis=0)
# combined_molecular_data.to_csv('data/downloads/cbioportal/rna_studies/rna_molecular_data.csv', index=False)

# combined_molecular_data = pd.read_csv('data/downloads/cbioportal/rna_studies/rna_molecular_data.csv')
# # Assertion:
# # Check if there is NA values
# combined_molecular_data.dropna(axis=1)