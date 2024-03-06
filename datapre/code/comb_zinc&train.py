import pandas as pd

# Load the train_set.csv file
xml_smiles_path = 'train_set.csv'
xml_smiles_df = pd.read_csv(xml_smiles_path)

# Load the existing zinc250k&non_spectra.csv file
zinc250k_path = 'zinc250k&non_spectra.csv'
zinc250k_df = pd.read_csv(zinc250k_path)

# Normalize column names to lowercase
xml_smiles_df.columns = [col.lower() for col in xml_smiles_df.columns]

# Remove the INCHI_KEY and INCHI_KEY_PREFIX columns
xml_smiles_df.drop(['inchi_key', 'inchi_key_prefix'], axis=1, inplace=True)

# Add default values for logP, qed, and SAS in xml_smiles_df
xml_smiles_df['logP'] = 0.0
xml_smiles_df['qed'] = 0.0
xml_smiles_df['SAS'] = 0.0

# Merge the two dataframes, ensuring that all entries from xml_smiles_df are included
merged_df = pd.merge(zinc250k_df, xml_smiles_df, on='smiles', how='outer', suffixes=('', '_y'))

# Drop the duplicate columns from the merge
merged_df.drop(columns=[col for col in merged_df if col.endswith('_y')], inplace=True)

# Fill any remaining NaN values with 0 (for new SMILES from xml_smiles_df that were not in zinc250k_df)
merged_df.fillna(0, inplace=True)

print(len(merged_df))

# Save the merged DataFrame to a new CSV file
merged_df.to_csv('nkn2024.csv', index=False)