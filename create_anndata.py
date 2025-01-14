import pandas as pd
import anndata as ad
import numpy as np

# This script should generate anndata files for a study
# Inputs
# 1. Molecular data -> adata.X
# 2. Meta data -> adata.obs
# 3. Output path
molecular_data_path = "data/downloads/cbioportal/tcga_pan_can_atlas/RNA/brca_tcga_pan_can_atlas_2018.csv"
meta_data_path = "data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv"
path_to_save_location = "data/inputs/PanCanAtlas_BRCA_transcriptome_RSEM.h5ad"
# Load the data
molecular_data = pd.read_csv(molecular_data_path)
meta_data = pd.read_csv(meta_data_path)

# Subset meta data to only include samples that are in molecular data
matched_meta_data = meta_data[meta_data['sampleId'].isin(molecular_data['sampleId'])]

# Check the number of unique samples
assert molecular_data['sampleId'].nunique() == matched_meta_data['sampleId'].nunique()

# Counts
counts = molecular_data.select_dtypes(include=np.number)
# Observations (in form of patient_id)
# Check if they are unique
assert molecular_data["patientId"].is_unique
observations = molecular_data["patientId"].values
# Features 
features = molecular_data.select_dtypes(include=np.number).columns

# Creating Anndata:
adata = ad.AnnData(counts)
adata.obs_names = observations
adata.var_names = features

# Add meta data to anndata object
columns_to_add = matched_meta_data.columns
for column in columns_to_add:
    adata.obs[column] = pd.Categorical(matched_meta_data[column])

# Make all variable names lower case
adata.obs.columns = map(str.lower, adata.obs.columns)
# Replace spaces with '_'
adata.obs = adata.obs.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)

# Rename some columns
column_mapping = {
    "cancertypeid": "cancer_type_id",
    "cancer_type_detailed": "cancer_type_detailed",
    "studyid": "study_id",
    "sampleid": "sample_id",
    "cancertype_name": "cancer_type_name",
    "consensus_ancestry": "genetic_ancestry",
    "pooled_consensus_ancestry": "genetic_ancestry",
    "consensus_ancestry": "genetic_ancestry_detailed",
    "race": "self_reported_ancestry",
    "sex": "sex"
}

adata.obs = adata.obs.rename(columns=column_mapping)

# Save
adata.write(path_to_save_location, compression="gzip")
