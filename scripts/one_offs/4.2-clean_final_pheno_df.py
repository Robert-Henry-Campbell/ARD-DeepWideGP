"""
Docstring for deepwidegp.ARD-DeepWideGP.scripts.one_offs.4.2-clean_final_pheno_df

goal: 
2. preserve only relevant columns for downstream analyses

"""

import pandas as pd
import os
os.chdir("/mnt/sdg/robert/deepwidegp/ARD-DeepWideGP/")

#load the df
final_pheno_df = pd.read_csv("data/4-integrate_MVP_UKB/2025_casecount/MVP_UKB_mapped_ards_MVP_GWAS_and_UKB_10000case.csv")

#select relevant columns
relevant_columns = [
'ICD10_category',
"ICD10_y",
"ICDO10_Phecode_mapped_flag",
"ICD10_x",
"ICD9",
"cause_level_1",
"cause_level_2",
"cause_level_3",
"cause_level_4",
'UKB_case_count',
"UKB_phenotype",
"mvp_match_flag",
"MVP_Trait",
"MVP_SubCategory",
"MVP_Description",
"MVP_num_samples.META",
"MVP_num_cases.META",
"MVP_num_controls.META",
"MVP_num_samples.AFR",
"MVP_num_cases.AFR",
"MVP_num_controls.AFR",
"MVP_num_samples.AMR",
"MVP_num_cases.AMR",
"MVP_num_controls.AMR",
"MVP_num_samples.EAS",
"MVP_num_cases.EAS",
"MVP_num_controls.EAS",
"MVP_num_samples.EUR",
"MVP_num_cases.EUR",
"MVP_num_controls.EUR",
]
final_pheno_df = final_pheno_df[relevant_columns]

#save the cleaned df
final_pheno_df.to_csv("data/4-integrate_MVP_UKB/2025_casecount/MVP_UKB_final_cleaned_pheno_df.csv", index=False)


pass