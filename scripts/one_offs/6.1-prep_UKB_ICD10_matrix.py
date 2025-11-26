import pandas as pd

#load the raw csv and the dictionary
map_df = pd.read_csv("ARD-DeepWideGP/data/3-identify UKB case count for ARDs/dna_nexus_wrangling/participant_data_dictionary.csv")
df = pd.read_csv("ARD-DeepWideGP/data/3-identify UKB case count for ARDs/dna_nexus_wrangling/first_occurance_24.11.25.csv")

#rename the columns 

# build mapping dict and apply rename
mapping = dict(zip(map_df['name'], map_df['title']))
df = df.rename(columns=mapping)

#create a pseudo identifier column
df['pseudo_id'] = range(1, len(df) + 1)


#save the new dataframe
df.to_csv("ARD-DeepWideGP/data/3-identify UKB case count for ARDs/dna_nexus_wrangling/ukb_icd10_matrix_2025_cleaned.csv", index=False)
pass
