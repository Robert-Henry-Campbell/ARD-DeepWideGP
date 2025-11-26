"""
Docstring for deepwidegp.ARD-DeepWideGP.scripts.one_offs.4.2-clean_final_pheno_df

goal: 
1. drop duplicate ICD10 codes from final MVP-UKB merged pheno df
2. preserve only relevant columns for downstream analyses

"""

import pandas as pd
import os
os.chdir("/mnt/sdg/robert/deepwidegp/ARD-DeepWideGP/")

#load the df
final_pheno_df = pd.read_csv("data/4-integrate_MVP_UKB/2025_casecount/MVP_UKB_mapped_ards_MVP_GWAS_and_UKB_10000case.csv")
#print the shape before
print("Shape before dropping duplicates:", final_pheno_df.shape)

#drop duplicate ICD10 codes
final_pheno_df = final_pheno_df.drop_duplicates(subset=['ICD10_explo'])
#print the shape after
print("Shape after dropping duplicates:", final_pheno_df.shape)

INCOMPLETE 
pass