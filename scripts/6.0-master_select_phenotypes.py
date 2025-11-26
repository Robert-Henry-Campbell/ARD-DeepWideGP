"""
What type of object do I want my snp list to be stored in?


function: create_UKB_pheno_casecount
inputs: 
1. cleaned UKB ICD10 pheno df
outputs:
1. UKB ICD10 pheno case count df

function: merge_GBD_ARD_UKB_MVP_thresholding
description: filters for phenos with ARD status in GBD, greather than (threshold) case coutn in UKB and present in MVP
inputs:
1. UKB ICD10 pheno case count df
2. ARD_GBD_ICD10_phecode_mapped_with_MVP df (avaliable in 2)
3. 



function: MVP_discovery_snps
Inputs: 
1. ICD10 pheno case count df
2. case count threshold (int)
optional: list of ICD10 codes (replaces thresholding)
3. output path (str)

Outputs: 
A XXXtype object containing of XXX-type files 
containing the SNPs of interest to be used by the deepGP pipeline. 


"""
