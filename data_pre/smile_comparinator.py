import pandas as pd

# Load the datasets
xml_data = pd.read_csv('datapre/data/xml_smiles_canonicalized.csv')
qm9_data = pd.read_csv('datapre/data/qm9_smiles_canonicalized.csv')
zinc250k_data = pd.read_csv('datapre/data/zinc250k_smiles_canonicalized.csv')

print("Length of xml data: ", len(xml_data))
print("Length of qm9 data: ", len(qm9_data))
print("Length of zinc data: ", len(zinc250k_data))

# Assuming the column containing SMILES is named 'smiles' in all datasets
xml_smiles_set = set(xml_data['SMILES'])
qm9_smiles_set = set(qm9_data['SMILES'])
zinc250k_smiles_set = set(zinc250k_data['SMILES'])

# Find SMILES in xml that are not in qm9 and not in zinc250k
unique_in_xml = xml_smiles_set - (qm9_smiles_set.union(zinc250k_smiles_set))

# Convert the set of unique SMILES back to a list and create a DataFrame
unique_smiles_df = pd.DataFrame(list(unique_in_xml), columns=['smiles'])

# Save the unique SMILES strings to a new CSV file
unique_smiles_df.to_csv('unique_xml_smiles.csv', index=False)

print(f"Found {len(unique_in_xml)} unique SMILES strings in the XML dataset.")
print("Unique SMILES strings saved to 'unique_xml_smiles.csv'")