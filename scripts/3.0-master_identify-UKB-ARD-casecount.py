"""
Docstring for deepwidegp.ARD-DeepWideGP.scripts.3.0-master_identify-UKB-ARD-casecount
run 3.1 
python scripts/3.1-create_UKB_casecount.py \
  --input data/3-identify-UKB-casecount-ARDs/dna_nexus_wrangling/ukb_icd10_matrix_2025_cleaned.csv \
  --output-dir data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025 \
  --output-name UKB_2025_casecount.csv

  
run 3.2 
python scripts/3.2-merge_GBD_ARD_UKB.py \
  --casecounts data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025/UKB_2025_casecount.csv \
  --ard-icd10 data/3-identify-UKB-casecount-ARDs/bothsex_ARD_GBD_ICD10_exploded_20250911_201039.csv \
  --output-dir data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025 \
  --output-name GBD_ARD_UKB_2025_ICD0_casecounts.csv

"""