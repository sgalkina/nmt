from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
import numpy as np
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')  
import matplotlib.pyplot as plt
from IPython.display import Image

### PICK KIND OF DATASAT
dataset_fil = 'valid'

vae_recons = pd.read_csv(f'{dataset_fil}_smiles_rec.csv')
nf_recons = pd.read_csv(f'OLD/{dataset_fil}_smiles_rec.csv')

print(vae_recons.shape)
print(nf_recons.shape)

original_smiles_column = 'Original SMILES'
reconstructed_smiles_column = 'Reconstructed SMILES' 

orig_smiles = vae_recons[original_smiles_column]
recon_smiles = vae_recons[reconstructed_smiles_column]
recon_NF_smiles = nf_recons[reconstructed_smiles_column]


molecule_to_visulize = 318

orig_smiles_rec = orig_smiles[molecule_to_visulize]
smi_nf = recon_NF_smiles[molecule_to_visulize]
smi = recon_smiles[molecule_to_visulize]

m = Chem.MolFromSmiles(smi,sanitize=False)
if m is None:
  print('invalid')
img = Draw.MolToImage(m)

fragments = smi.split('.')
molecules = []

for fragment in fragments:
    mol = Chem.MolFromSmiles(fragment, sanitize=False)
    if mol is not None:
        molecules.append(mol)
    else:
        print(f'Invalid fragment: {fragment}')

# Find the largest fragment by the number of heavy atoms
if molecules:
    largest_mol = max(molecules, key=lambda mol: mol.GetNumHeavyAtoms())

    # Visualize the largest fragment
    img = Draw.MolToImage(largest_mol)
    img.show()
else:
    print('No valid molecules to visualize.')

mol = Chem.MolFromSmiles(orig_smiles_rec)


if mol is not None:
    Chem.AllChem.Compute2DCoords(mol)
    
    # Draw the molecule using RDKit
    drawer = Draw.MolDraw2DCairo(300, 300)  # Drawing area size in pixels
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Convert the drawing to a PNG image
    png_data = drawer.GetDrawingText()
    
    # Display the image using Matplotlib
    plt.figure(figsize=(5, 5))
    with open("mol.png", "wb") as f:
        f.write(png_data)
    plt.imshow(plt.imread("mol.png"))
    plt.axis('off')  # Hide the axes
    plt.show()
else:
    print("Invalid SMILES string.")



### FOR FINDING MEAN AND STD'S
def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:  
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    else:
        return None

# Compute similarities
similarities = []
for index, row in vae_recons.iterrows():
    orig_smiles = row[original_smiles_column]
    recon_smiles = row[reconstructed_smiles_column]
    
    orig_fp = get_fingerprint(orig_smiles)
    recon_fp = get_fingerprint(recon_smiles)
    
    if orig_fp is not None and recon_fp is not None:
        similarity = DataStructs.TanimotoSimilarity(orig_fp, recon_fp)
        similarities.append(similarity)

# Compute mean and standard deviation
if similarities:
    mean_similarity = np.mean(similarities)
    std_similarity = np.std(similarities)
    print("Mean Similarity:", mean_similarity)
    print("Standard Deviation:", std_similarity)
else:
    print("No valid similarities computed.")




"""
index_to_view = 2
hmdb['Match'] = hmdb[original_smiles_column] == hmdb[reconstructed_smiles_column]
matching_count = hmdb['Match'].sum()
print(f"Number of identical SMILES pairs: {matching_count}")
print(f"Percentage of match: {matching_count / len(hmdb) * 100:.2f}%")

dataset_fil = 'train'
hmdb = pd.read_csv(f'{dataset_fil}_smiles_rec.csv')
index_to_view = 2

def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:  # Check if molecule conversion is successful
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    else:
        return None

# Compute similarities
similarities = []
for orig, recon in zip(original_smiles_column, reconstructed_smiles_column):
    orig_fp = get_fingerprint(orig)
    recon_fp = get_fingerprint(recon)
    if orig_fp is not None and recon_fp is not None:
        similarity = DataStructs.TanimotoSimilarity(orig_fp, recon_fp)
        similarities.append(similarity)

# Compute mean and standard deviation
if similarities:
    mean_similarity = np.mean(similarities)
    std_similarity = np.std(similarities)
    print("Mean Similarity:", mean_similarity)
    print("Standard Deviation:", std_similarity)
else:
    print("No valid similarities computed.")

if index_to_view < len(hmdb):
    original_smiles = hmdb.loc[index_to_view, original_smiles_column]
    reconstructed_smiles = hmdb.loc[index_to_view, reconstructed_smiles_column]

    # Create RDKit molecule objects from SMILES strings
    original_molecule = Chem.MolFromSmiles(original_smiles)
    reconstructed_molecule = Chem.MolFromSmiles(reconstructed_smiles)

    # Check if molecules are successfully created
    if original_molecule and reconstructed_molecule:
        combined_molecules = [original_molecule, reconstructed_molecule]
        legends = [f'Orig: {original_smiles}', f'Recon: {reconstructed_smiles}']

        # Visualize the molecules
        img = Draw.MolsToGridImage(combined_molecules, molsPerRow=2, subImgSize=(200, 200), legends=legends)
        img.show()

        # Save the combined image
        img.save(f'Reconstructions_Visualize/{dataset_fil}/combined_molecules_{index_to_view+1}.png')
        print(f'Saved reconstruction image as combined_molecules_{index_to_view+1}.png')
    else:
        print("Error: Failed to create molecules from SMILES.")
else:
    print(f"Error: Index {index_to_view} is out of bounds. Maximum index is {len(hmdb)-1}.") """