import pandas as pd
import anndata as ad
import numpy as np

# This script should generate anndata files for a study
# Inputs
# 1. Molecular data -> adata.X
# 2. Meta data -> adata.obs
# 3. Output path
molecular_data_path = "data/downloads/cbioportal/tcga_firehose/methylation/brca_tcga.csv"

meta_data_path = "data/downloads/cbioportal/tcga_pan_can_atlas/meta_tcga_pan_can_atlas_expression.csv"

path_to_save_location = "data/inputs/PanCanAtlas_BRCA_raw_BETA_subtypeNAremoved.h5ad"
# Load the data
molecular_data = pd.read_csv(molecular_data_path)
meta_data = pd.read_csv(meta_data_path)
# meta_data[meta_data["studyId"] == "brca_tcga_pan_can_atlas_2018"]["CANCER_TYPE_DETAILED"].unique()

# Subset meta data 
matched_meta_data = meta_data[meta_data["sampleId"].isin(molecular_data["sampleId"])]

# Subset molecular data
matched_molecular_data = molecular_data[molecular_data["sampleId"].isin(matched_meta_data["sampleId"])]

# Check the number of unique patients
assert matched_molecular_data["patientId"].nunique() == matched_meta_data["patientId"].nunique()
assert matched_molecular_data["sampleId"].nunique() == matched_meta_data["sampleId"].nunique()

# Counts
counts = matched_molecular_data.select_dtypes(include=np.number)
# Observations (in form of patient_id)
# Check if they are unique
assert matched_meta_data["patientId"].is_unique
observations = matched_meta_data["patientId"].values
# Drop so its not doubled
matched_meta_data = matched_meta_data.drop(columns="patientId")
# Features 
features = matched_molecular_data.select_dtypes(include=np.number).columns

# Creating Anndata:
adata = ad.AnnData(counts)
adata.obs_names = observations
adata.var_names = features

# Remove Na values
columns_with_na = adata.var_names[np.isnan(adata.X).any(axis=0)]
adata = adata[:, ~adata.var_names.isin(columns_with_na)]
# Remove negative values (Set them to zero)
# adata.X[adata.X < 0] = 0

# Add meta data to anndata object
columns_to_add = matched_meta_data.columns
for column in columns_to_add:
    adata.obs[column] = pd.Categorical(matched_meta_data[column])

# Convert age to numeric
adata.obs["AGE"] = adata.obs["AGE"].cat.codes.astype(int)

# Make all variable names lower case
adata.obs.columns = map(str.lower, adata.obs.columns)
# Replace spaces with '_'
adata.obs = adata.obs.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)
adata.obs["cancer_type_detailed"].unique()

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

# BRCA
# Remove all NAs
adata = adata[adata.obs.dropna(subset="subtype").index]

# Include subtype_pooled
adata.obs["subtype_pooled"] = adata.obs["subtype"].apply(
    lambda x: "Basal" if x == "BRCA_Basal" else "non-Basal"
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
