# This script needs a working conda environment with pybioportal installed
# Pybioportal
import pybioportal 
from pybioportal import studies as stdy
from pybioportal import sample_lists as sl
from pybioportal import clinical_data as cd
# Standard libraries
import pandas as pd
import numpy as np

# Get all TCGA studies 
# 1. Get all studies with keyword: TCGA
# 2. Save the meta data of the studies
tcga_studies = stdy.get_all_studies(keyword="TCGA", projection='DETAILED') 
tcga_rna_studies = tcga_studies[tcga_studies['mrnaRnaSeqV2SampleCount'] != 0] 
# tcga_rna_studies.to_csv('data/downloads/cbioportal/tcga_studies_with_rna.csv', index=False)

# Interesting studies are either from Firehose Legacy or PanCancerAtlas
# 1. Seperate those entries based on the origin
# 2. Keep working with the studies from the PanCancerAtlas
# (PanCancer studies dont have any other data than RNA)
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
    
    # study_id = tcga_rna_studies_pan['studyId'].values[0]
    # Filter important information
    attributes = ['name', 'studyId', 'cancerTypeId', 'cancerType_name'] # cancerTypeId is a short name for the cancer
    study_information = pybioportal.studies.fetch_studies(study_ids=[study_id], projection='DETAILED')
    study_information = study_information[attributes]

    available_sample_lists = sl.get_all_sample_lists_in_study(study_id=study_id, projection='DETAILED')
    # available_sample_lists.sampleListId
    rna_sample_list = available_sample_lists[available_sample_lists['category'] == 'all_cases_with_mrna_rnaseq_data']
    # Get the list of patients in the column 'sampleIds'
    rna_samples = rna_sample_list['sampleIds'][0]

    # Assertion: 
    # Check if 'sampleCount' is not empty
    assert len(rna_samples) > 0, \
        f'StudyId: {study_id} has no samples available.'

    # Check if there is dublicated samples with '-01' extension
    assert len(rna_samples) == len(set(rna_samples)), \
        f'StudyId: {study_id} has dublicated samples.'
    
    # Without extension
    # Remove the last segment (e.g., "-01", "-02", etc.)
    stripped_samples = [sample.rsplit('-', 1)[0] for sample in rna_samples]
    assert len(rna_samples) == len(stripped_samples), \
        f'After removing extension samples are not unique'

    # Combine information to dataframe 
    # 1. Add 'study_id' column to merge information 
    # 2. Add 'study_information' (merge with 'meta_df' on 'studyId')
    meta_df = pd.DataFrame(data=rna_samples, columns=['sampleId'])
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
            'ETHNICITY', 'SUBTYPE']
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

    # Check for dubplicated samples
    if not meta_df.sampleId.is_unique:
        print(f'StudyId: {study_id} has replicated samples.')
    
    assert meta_df.sampleId.is_unique, \
        f'StudyId: {study_id} has replicated samples.'

    # Combine all studies
    meta_all = pd.concat([meta_all, meta_df])

# meta_all.to_csv('data/downloads/cbioportal/tcga_studies_with_rna.csv', index=False)

# Done: 
# Checked if samples are unique per study
# Checked if patients are unique per study

# Further quality check 
# 1. Check if samples are unique across studies 
# 2. Check if patients are unnique across studies (possibility that one patients has more than one cancer)
# meta_all = pd.read_csv('data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv')

assert meta_all['sampleId'].is_unique, \
    'Samples are not unique across studies.'

assert meta_all['patientId'].is_unique, \
    'There is patients across studies.'

# Add genetic ancestry (genetically inferred by doi.org/10.1016/j.ccell.2020.04.012)
# Excel file with the ancestry
xls = pd.ExcelFile('../data/downloads/cbioportal/1-s2.0-S1535610820302117-mmc2.xlsx') 
# Sheet with ancestry                               
ancestry_sheet = pd.read_excel(xls, "S1 Calls per Patient", header=1)            
ancestry_sheet = ancestry_sheet[['patient', 'tumor_type', 'consensus_ancestry']]
ancestry_sheet.rename(columns={'patient': 'patientId'}, inplace=True)

# Merge 'ancestry_sheet' with 'meta' dataframe
# This procedure loses patients without genetic ancestry
meta_all = meta_all.merge(ancestry_sheet, on='patientId')

# Pooling admix into one ancestry
meta_all['pooled_consensus_ancestry'] = np.where(meta_all['consensus_ancestry'].str.contains('admix'), 
                                                 'admix',
                                                 meta_all['consensus_ancestry'])

# Add column to notify that RNA Seq V2 data is available
meta_all['RNASeqV2'] = True

# Do reformating
# 1. Make all colnames lower case
# meta_all.columns = map(str.lower, meta_all.columns)

# Save meta dataframe
meta_all.to_csv('../data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv', index=False)

# meta = pd.read_csv('data/downloads/cbioportal/rna_studies/tcga_pan_studies_with_rna.csv')




