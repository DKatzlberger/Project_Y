# This script needs a working conda environment with pybioportal installed
import pybioportal 
from pybioportal import studies as stdy
from pybioportal import sample_lists as sl
from pybioportal import clinical_data as cd
import pandas as pd

# Get all TCGA studies 
# 1. Get all studies with keyword: TCGA
# 2. Get studies with RNA Seq V2
# 3. Save the meta data of the studies
tcga_studies = stdy.get_all_studies(keyword="TCGA", projection='DETAILED') 
tcga_rna_studies = tcga_studies[tcga_studies['mrnaRnaSeqV2SampleCount'] != 0] 
tcga_rna_studies.to_csv('data/downloads/cbioportal/tcga_studies_with_rna.csv', index=False)

# Interesting studies are either from Firehose Legacy or PanCancerAtlas
# 1. Seperate those entries based on the origin
# 2. Keep working with the studies from the PanCancerAtlas
tcga_rna_studies_pan = tcga_rna_studies[tcga_rna_studies['name'].str.contains('TCGA, PanCancer Atlas')]
tcga_rna_studies_fire = tcga_rna_studies[tcga_rna_studies['name'].str.contains('TCGA, Firehose Legacy')]

# Retrieve  patients in the PanCancerAtlas
# Loop over all studies in 'tcga_rna_studies_pan'
# 1. Get information about the study
# 2. Retrieve all patient lists
# 3. Filter for RNA Seq V2 ('category' == 'all_cases_with_mrna_rnaseq_data')
# 4. Get patients with RNA Seq V2 data
meta_all = pd.DataFrame()
for study_id in tcga_rna_studies_pan['studyId']:
    print(study_id)
    
    # study_id = tcga_rna_studies_pan['studyId'][1]
    # Filter important information
    attributes = ['name', 'studyId', 'cancerTypeId', 'cancerType_name'] # cancerTypeId is a short name for the cancer
    study_information = pybioportal.studies.fetch_studies(study_ids=[study_id], projection='DETAILED')
    study_information = study_information[attributes]

    available_sample_lists = sl.get_all_sample_lists_in_study(study_id=study_id, projection='DETAILED')
    rna_patient_list = available_sample_lists[available_sample_lists['category'] == 'all_cases_with_mrna_rnaseq_data']
    # Get the list of patients in the column 'sampleIds'
    rna_patients = rna_patient_list['sampleIds'][0]

    # Assertion: 
    # Check if 'sampleCount' is not empty
    assert len(rna_patients) > 0, \
        f'StudyId: {study_id} has no patients available.'

    # Check if there is dublicates
    assert len(rna_patients) == len(set(rna_patients)), \
        f'StudyId: {study_id} has dublicated patients.'

    # Combine information to dataframe 
    # 1. Add 'study_id' column to merge information 
    # 2. Add 'study_information' (merge with 'meta_df' on 'studyId')
    meta_df = pd.DataFrame(data=rna_patients, columns=['sampleId'])
    meta_df['studyId'] = study_id
    meta_df = meta_df.merge(study_information, on='studyId')

    # Get all the sample info
    # 1. Filter important information
    # 2. Add 'sample_information' (merge with 'meta_df' on 'sampleId')
    attributes = ['sampleId', 'patientId', 'CANCER_TYPE_DETAILED']
    sample_information = cd.fetch_all_clinical_data_in_study(study_id=study_id, clinical_data_type='SAMPLE')
    sample_information = sample_information[attributes]
    # sample_information.columns.name = None

    meta_df = meta_df.merge(sample_information, on='sampleId')

    # Get all the patient info 
    # 1. Filter important information
    # 2. Add 'patient_information' (merge with 'meta_df' on 'patientId')
    attributes = ['patientId', 'AGE', 'SEX', 'RACE', \
            'ETHNICITY', 'GENETIC_ANCESTRY_LABEL', 'SUBTYPE']
    patient_information = cd.fetch_all_clinical_data_in_study(study_id=study_id, clinical_data_type='PATIENT')  
    # Insert loop to try to get each item in 'attribute'
    # 1. Try to retrieve attribute
    # 2. Except print statement there is not attribute
    patient_information_all = pd.DataFrame()
    for attribute in attributes:
        try:
            attribute_information = patient_information[attribute]
            patient_information_all = pd.concat([patient_information_all, attribute_information], axis=1)
        except:
            print(f'StudyId: {study_id} has no patient information: {attribute}.')

    # Merge
    meta_df = meta_df.merge(patient_information_all, on='patientId')

    # Some studies have replicates in the data
    if not meta_df.patientId.is_unique:
        print(f'StudyId: {study_id} has replicated patients, removing them.')

    # Remove dublicated patients
    # 1. Find duplicated patients
    # 2. Remove duplicated patient
    # 3. Check if removed
    duplicates = meta_df[meta_df.duplicated('patientId')].index
    meta_df = meta_df[~meta_df.index.isin(duplicates)]

    assert meta_df.patientId.is_unique, \
        f'StudyId: {study_id} has replicated patients.'

    if not meta_df.sampleId.is_unique:
        print(f'StudyId: {study_id} has replicated samples.')
        
    # assert meta_df.sampleId.is_unique, \
    #     f'StudyId: {study_id} has replicated samples.'

    # Do reformating
    # 1. Make all colnames lower case
    meta_df.columns = map(str.lower, meta_df.columns)
    # meta_df.to_csv(f'data/downloads/cbioportal/{study_id}.csv', index=False)

    # Combine all studies
    meta_all = pd.concat([meta_all, meta_df])

meta_all.to_csv('data/downloads/cbioportal/tcga_studies_with_rna.csv', index=False)

meta = pd.read_csv('data/downloads/cbioportal/tcga_studies_with_rna.csv')

meta['sampleid'].is_unique
