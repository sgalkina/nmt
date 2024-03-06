from rdkit import Chem
import pandas as pd
from sklearn.model_selection import train_test_split

# Load the xml_smiles_unique dataset
xml_smiles_unique_path = 'xml_smiles_unique.csv'
xml_smiles_unique_df = pd.read_csv(xml_smiles_unique_path)

# Load the SDF file
sdf_path = 'structures.sdf'
sdf_supplier = Chem.SDMolSupplier(sdf_path)

# Create a DataFrame to store SMILES and INCHI_KEY
sdf_data = []
for mol in sdf_supplier:
    if mol is not None:
        sdf_data.append({
            'SMILES': Chem.MolToSmiles(mol),
            'INCHI_KEY': mol.GetProp('INCHI_KEY') if mol.HasProp('INCHI_KEY') else None
        })

sdf_df = pd.DataFrame(sdf_data)

# Merge to get INCHI_KEY for each SMILES in xml_smiles_unique_df
merged_df = pd.merge(xml_smiles_unique_df, sdf_df, left_on='smiles', right_on='SMILES', how='left')

# After merging, you may want to drop the redundant 'SMILES' column from sdf_df
merged_df.drop(columns=['SMILES'], inplace=True)

# Extract the first 14 characters of the INCHI_KEY to ensure unique chemical identity
merged_df['INCHI_KEY_PREFIX'] = merged_df['INCHI_KEY'].apply(lambda x: x[:14] if pd.notnull(x) else None)

# Group by INCHI_KEY_PREFIX to ensure same molecules are not split across sets
grouped = merged_df.groupby('INCHI_KEY_PREFIX')

# Splitting groups into train+valid and test sets
train_valid_groups, test_groups = train_test_split(list(grouped), test_size=0.2, random_state=42)

# Convert groups back to DataFrames
train_valid_df = pd.concat([group[1] for group in train_valid_groups], ignore_index=True)
test_df = pd.concat([group[1] for group in test_groups], ignore_index=True)

# Further split train_valid_df into train and validation sets
train_df, valid_df = train_test_split(train_valid_df, test_size=0.25, random_state=42)  

# Save the splits to CSV files
train_df.to_csv('train_set.csv', index=False)
valid_df.to_csv('valid_set.csv', index=False)
test_df.to_csv('test_set.csv', index=False)

print(f"Train set size: {len(train_df)}")
print(f"Validation set size: {len(valid_df)}")
print(f"Test set size: {len(test_df)}")