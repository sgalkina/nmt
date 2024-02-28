import pandas as pd
from rdkit import Chem

data_name = 'zinc250k'
# data_name = 'xml_smiles'
# data_name= 'qm9_smiles'

data_path = f'{data_name}_smiles.csv'
data = pd.read_csv(data_path)

### For canonilization
def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, canonical=True)
    else:
        return None

# Apply the canonicalization function to the 'smiles' column
data['SMILES'] = data['SMILES'].apply(canonicalize_smiles)

# Save the result to a new CSV file
data.to_csv(f'{data_name}_smiles_canonicalized.csv', index=False)

print("Canonicalized SMILES strings saved to 'qm9_smiles_canonicalized.csv'")