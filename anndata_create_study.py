import pandas as pd
import anndata as ad
import numpy as np
import os
# This script should generate anndata files for a study
# Inputs
# 1. Molecular data -> adata.X
# 2. Meta data      -> adata.obs
molecular_data_path = "data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients/Molecular_data.csv"
meta_data_path      = "data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients/Meta_data.csv"
# 3. Output path
path_to_save_location = "data/downloads/GEO/Population-specific Mutation Patterns in Breast Tumors from African American, European American, and Kenyan Patients"

# Load the data
molec = pd.read_csv(molecular_data_path)
meta  = pd.read_csv(meta_data_path)

# Match data
sample_name   = "sample_title"
matched_meta  = meta[meta[sample_name].isin(molec[sample_name])]
matched_molec = molec[molec[sample_name].isin(matched_meta[sample_name])]
# Check the number of unique patients
assert matched_molec[sample_name].nunique() == matched_meta[sample_name].nunique()

# Counts
counts = matched_molec.select_dtypes(include=np.number)
# Change negative values to 0
counts = counts.clip(lower=0)
# Remove all NA columns
NA_columns = counts.columns[counts.isna().any(axis=0)].tolist()
counts     = counts.drop(columns=NA_columns)
print(f"Dropping {len(NA_columns)} columns because they contain NA values.")
# Remove all 0 columns
zero_columms = counts.columns[(counts == 0).all(axis=0)]
counts       = counts.drop(columns=zero_columms)
print(f"Dropping {len(zero_columms)} columns because they have no variance.")

# Observations 
obs = matched_meta[sample_name].values
# Check uniquness
assert len(obs) == len(set(obs)), "Error: Duplicated observations!"

# Features 
features = counts.columns.values
# Check uniquness
assert len(features) == len(set(features)), "Error: Duplicated features!"

# AnnData
adata = ad.AnnData(counts)
adata.obs_names = obs
adata.var_names = features

# Add meta columns
columns_to_add = matched_meta.drop(columns=sample_name).columns
for column in columns_to_add:
    adata.obs[column] = pd.Categorical(matched_meta[column])

# Make all variable names lower case
adata.obs.columns = map(str.lower, adata.obs.columns)
# Replace spaces with '_'
adata.obs = adata.obs.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)

# Save
save_name = os.path.join(path_to_save_location, "GSE225846_raw_RSEM.h5ad")
adata.write(save_name, compression="gzip")


# Convert age to numeric
adata.obs["AGE"] = adata.obs["AGE"].cat.codes.astype(int)

# Rename some columns
column_mapping = {
    "cancertypeid": "cancer_type_id",
    "cancer_type_detailed": "cancer_type_detailed",
    "studyid": "study_id",
    "sampleid": "sample_id",
    "patientid": "patient_id",
    "cancertype_name": "cancer_type_name",
    "consensus_ancestry": "genetic_ancestry",
    "pooled_consensus_ancestry": "genetic_ancestry",
    "consensus_ancestry": "genetic_ancestry_detailed",
    "race": "self_reported_ancestry",
    "sex": "sex"
}

adata.obs = adata.obs.rename(columns=column_mapping)

# Custom modifications

# Remove all NAs in subtype
adata = adata[adata.obs.dropna(subset="subtype").index]

# BRCA
# Include subtype_pooled
# adata.obs["subtype_pooled"] = adata.obs["subtype"].apply(
#     lambda x: "Basal" if x == "BRCA_Basal" else "non-Basal"
# )

# UCEC
adata.obs["subtype"] = adata.obs["subtype"].apply(
    lambda x: "CN-high" if x == "UCEC_CN_HIGH" else "non-CN-high"
)

# UCEC
# Remove parts of hte name after /
# adata.obs["cancer_type_detailed"] = adata.obs["cancer_type_detailed"].str.split("/").str[0]
# adata = adata[adata.obs.dropna(subset="subtype").index]
# # # Include subtype_pooled
# adata.obs["subtype_pooled"] = adata.obs["subtype"].apply(
#     lambda x: "cn_high" if x == "UCEC_CN_HIGH" else "non_cn_high"
# )

# Save
adata.write(path_to_save_location, compression="gzip")
