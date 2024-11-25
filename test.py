import pandas as pd
admix = pd.read_csv('data/combined_runs/TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_EUR_to_ADMIX/Weights.csv')
amr = pd.read_csv('data/combined_runs/TCGA_Breast_Invasive_Ductal_Carcinoma_vs_Breast_Invasive_Lobular_Carcinoma_EUR_to_AFR/Weights.csv')

admix.columns.isin(amr.columns).all()