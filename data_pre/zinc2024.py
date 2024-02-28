import pandas as pd

# Load the datasets
zinc250k_df = pd.read_csv('zinc250k.csv')  
xml_smiles_df = pd.read_csv('unique_xml_smiles.csv')  

# Rename the column in xml_smiles_df to match zinc250k_df, if needed
# xml_smiles_df.columns = ['smiles']  # Uncomment if the column name in xml_smiles_df is different

# Check for smiles in xml_smiles_df that are not in zinc250k_df
new_smiles = ~xml_smiles_df['smiles'].isin(zinc250k_df['smiles'])

# Create a DataFrame for the new smiles with 0 for logP, qed, and SAS
new_smiles_df = xml_smiles_df[new_smiles].copy()
new_smiles_df['logP'] = 0
new_smiles_df['qed'] = 0
new_smiles_df['SAS'] = 0

# Combine the two DataFrames
combined_df = pd.concat([zinc250k_df, new_smiles_df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('data2024.csv', index=False)

print(len(xml_smiles_df))
print(len(zinc250k_df))
print(len(combined_df))
print("Combined dataset saved as combined_dataset.csv")