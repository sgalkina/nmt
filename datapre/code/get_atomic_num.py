from rdkit import Chem
import pandas as pd

# Path to your CSV file
data = '/Users/nikolaikjaernielsen/Desktop/nkn2024.csv'

# Load the dataset from the CSV file
df = pd.read_csv(data)

unique_atomic_numbers = set()

for index, row in df.iterrows():
    molecule = Chem.MolFromSmiles(row['smiles'])
    if molecule is not None:
        for atom in molecule.GetAtoms():
            unique_atomic_numbers.add(atom.GetAtomicNum())

unique_atomic_numbers = sorted(list(unique_atomic_numbers))

print("Unique atomic numbers in the dataset:", unique_atomic_numbers)

"""
nkn2024.csv unique atomic numbers: 
Unique atomic numbers in the dataset: [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 
18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 
42, 44, 46, 47, 48, 50, 51, 52, 53, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 71, 
72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 90, 96]
"""