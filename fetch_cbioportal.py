# This script needs a working conda environment with pybioportal installed
import pybioportal 
import pandas as pd



# Study id
study_id = "brca_tcga"

# Clinical data (meta data)
attributes = ['patientId', 
              'studyId', 
              'AGE', 
              'SEX', 
              'RACE', 
              'SAMPLE_COUNT', 
              'ETHNICITY'
              ]
# Fetch clinical data for study id (=Meta data)
patient_data = pybioportal.clinical_data.fetch_all_clinical_data_in_study(study_id = study_id, 
                                                                           clinical_data_type="PATIENT",
                                                                           # attribute_ids=attributes
                                                                           )        
# Fetch sample data
sample_data = pybioportal.clinical_data.fetch_all_clinical_data_in_study(study_id = study_id, 
                                                                           clinical_data_type="SAMPLE",
                                                                           # attribute_ids=attributes
                                                                           )


# Remove column name
patient_data.columns.name = None
sample_data.columns.name = None

# Rename column names


# Check uniquness of patients
assert patient_data['patientId'].nunique() == len(patient_data), \
    f'Patients not unique in study: {study_id}'


# Check if patients are not in any other TCGA studies
# TODO - check if patients are unique in the cancer study (even though other cancers are not in the same tissue)

# Merge patient and sample data
meta = patient_data.merge(sample_data, on='patientId')

# Check how many patients have replicates 
meta['SAMPLE_COUNT'] = meta['SAMPLE_COUNT'].astype(int)
replicated_patient = meta[meta['SAMPLE_COUNT'] > 1]

# Filter columns by selected attributes
meta_filtered = meta[attributes]

# 


