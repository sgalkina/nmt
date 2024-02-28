import pandas as pd

# Load the datasets
qm9_data = pd.read_csv('qm9_smiles_canonicalized.csv')
zinc250k_data = pd.read_csv('zinc250k_smiles_canonicalized.csv')
xml_data = pd.read_csv('xml_smiles_canonicalized.csv')

# Concatenate the datasets
combined_data = pd.concat([qm9_data, zinc250k_data, xml_data], ignore_index=True)

# Optional: Remove duplicate SMILES strings, if desired
combined_data.drop_duplicates(subset='SMILES', inplace=True)

# Save the combined dataset to a new CSV file
combined_data.to_csv('combined_smiles_dataset.csv', index=False)

print("Combined dataset saved to 'combined_smiles_dataset.csv'. Total entries:", len(combined_data))

### USE THIS TO LOAD COMBINED DATA
### SHOULD HAVE 385.501 ROWS
combined_dataset = pd.read_csv('combined_smiles_dataset.csv', low_memory=False)
print(len(qm9_data))