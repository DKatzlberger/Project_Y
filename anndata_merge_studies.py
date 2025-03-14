import pandas as pd
import anndata as ad
import numpy as np

study_I_path = "data/downloads/cbioportal/tcga_pan_can_atlas/expression/lusc_tcga_pan_can_atlas_2018.csv"
study_II_path= "data/downloads/cbioportal/tcga_pan_can_atlas/expression/luad_tcga_pan_can_atlas_2018.csv"
meta_data_path = "data/downloads/cbioportal/tcga_pan_can_atlas/meta_tcga_pan_can_atlas_expression.csv"
# Name
path_to_save_location = "data/inputs/PanCanAtlas_LUSC_LUAD_raw_RSEM_subtypeNAremoved.h5ad"
# Load the data
study_I = pd.read_csv(study_I_path)
study_II = pd.read_csv(study_II_path)
meta_data = pd.read_csv(meta_data_path)

# Concat rowwise
study_merged = pd.concat([study_I, study_II])

# Subset meta data to only include samples that are in molecular data
matched_meta_data = meta_data[meta_data['sampleId'].isin(study_merged['sampleId'])]

# Check the number of unique samples
assert study_merged['sampleId'].nunique() == matched_meta_data['sampleId'].nunique()

#Counts
counts = study_merged.select_dtypes(include=np.number)
# Observations (in form of patient_id)
# Check if they are unique
assert matched_meta_data["patientId"].is_unique
observations = matched_meta_data["patientId"].values
# Drop so its not doubled
matched_meta_data = matched_meta_data.drop(columns="patientId")
# Features 
features = study_merged.select_dtypes(include=np.number).columns

# Creating Anndata:
adata = ad.AnnData(counts)
adata.obs_names = observations
adata.var_names = features

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
# Remove all NAs
adata = adata[adata.obs.dropna(subset="subtype").index]

# # Include subtype_pooled
# adata.obs["subtype_pooled"] = adata.obs["subtype"].apply(
#     lambda x: "basal" if x == "BRCA_Basal" else "non_basal"
# )

# Save
adata.write(path_to_save_location, compression="gzip")
