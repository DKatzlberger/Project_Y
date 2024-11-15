import pandas as pd
import anndata as ad
import numpy as np

# Meta data with the original ancestries
meta_original = pd.read_csv('data/downloads/tcga_csv/tcga_studies_meta_original_ancestry.csv')
meta_original = meta_original.rename(columns={'consensus_ancestry': 'genetic_ancestry_detailed'})

# Load latest data of Master thesis project
meta = pd.read_csv('data/downloads/tcga_csv/tcga_studies_meta_intersection_adjusted_withnewcolumns.csv', index_col=[0])
data = pd.read_csv('data/downloads/tcga_csv/hm450_merged_betavalues.csv', index_col=[0])

# Add non-pooled ancestry (original ancestries) to meta
meta = meta.merge(meta_original[['sampleId', 'genetic_ancestry_detailed']], how='left', on='sampleId')

# Replace all spaces with underscore
meta = meta.replace(' ', '_', regex=True)

# Check if they are the same length
assert meta.shape[0] == data.shape[0], 'Dataframes do not have the same length'
# Check if they have the same samples
assert meta.sampleId.isin(data.sampleId).all(), 'Not the same samples'

# Bring them in the same order
meta = meta.sort_values(by='sampleId').reset_index(drop=True)
data = data.sort_values(by='sampleId').reset_index(drop=True)

# Check the order
assert meta.sampleId.equals(data.sampleId), 'Meta and data not in right order'

# Preparing for Anndata
# Gene counts of the data
counts = data.select_dtypes(include=np.number).values
# Observation names (sampleIds)
obs_names = data.sampleId.values
# Feature names (Gene names)
var_names = data.drop(columns='sampleId').columns

# Check dimensions
# Observations
assert counts.shape[0] == obs_names.shape[0], f'Counts {counts.shape[0]} and observations {obs_names.shape[0]}'
# Features
assert counts.shape[1] == var_names.shape[0], f'Counts {counts.shape[1]} and features {var_names.shape[0]}'

# Create anndata object
adata = ad.AnnData(counts)
adata.obs_names = obs_names
adata.var_names = var_names

# Add meta information
adata.obs['cancer_type'] = pd.Categorical(meta.cancerId)
adata.obs['cancer_type_detailed'] = pd.Categorical(meta.cancer_type_detailed)
adata.obs['study_id'] = pd.Categorical(meta.studyId)
adata.obs['genetic_ancestry'] = pd.Categorical(meta.consensus_ancestry)
adata.obs['genetic_ancestry_detailed'] = pd.Categorical(meta.genetic_ancestry_detailed)
adata.obs['self_reported_ancestry'] = pd.Categorical(meta.race)
adata.obs['sex'] = pd.Categorical(meta.sex)


# Save anndata
adata.write('TCGA_transcriptome_RSEM.h5ad', compression="gzip")