# Library
import pandas as pd

# Exchange meta data of patients of Pan Cancer analysis (has most information per patient)
pan_cancer_meta = pd.read_csv("data/downloads/cbioportal/tcga_pan_can_atlas/meta_tcga_pan_can_atlas_expression.csv")
pan_cancer_attributes = set(pan_cancer_meta.columns.values)

# Other project meta (e.g. Firehose)
fire_hose_meta = pd.read_csv("data/downloads/cbioportal/tcga_firehose/meta_tcga_firehose_protein.csv")
fire_hose_attributes = set(fire_hose_meta.columns.values)
# Missing attriubutes in other study
attribute_diff = pan_cancer_attributes - fire_hose_attributes

# Add missing attributes of the Pan Cancer by sample id
subset = list(attribute_diff)
subset.append("patientId")
to_add = pan_cancer_meta[subset]

# Merge with other project
merged_df = pd.merge(fire_hose_meta, to_add, on='patientId', how='inner')

# Remove rows with NA in attribute_diff
merged_df = merged_df.dropna(subset=list(attribute_diff))

# Save
merged_df.to_csv("data/downloads/cbioportal/tcga_firehose/meta_tcga_firehose_shared_protein.csv", index=False)
