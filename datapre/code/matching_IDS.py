import pandas as pd
from rdkit import Chem

# Load the datasets
structures_df = pd.read_csv('struc_databse_ids.csv')  # Adjust path as needed
xml_database_df = pd.read_csv('xml_database_ids.csv')  # Adjust path as needed

# Assuming 'database_ids' is the column name in both CSV files
# and 'smile' is the column name for the SMILES strings in structures.csv

# Find database_ids present in both dataframes
common_ids = pd.merge(structures_df, xml_database_df, on='DATABASE_ID')
# common_ids = common_ids.drop_duplicates(subset='DATABASE_ID')

print(len(structures_df))
print(len(xml_database_df))
print(len(common_ids))

### For canonilization
def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, canonical=True)
    else:
        return None

common_ids['SMILES'] = common_ids['SMILES'].apply(canonicalize_smiles)


# Save the result to a new CSV file
common_ids.to_csv('xml_smiles_canonicalized.csv', index=False)

print("Common database IDs with matching SMILES strings saved to 'xml_smiles.csv'")