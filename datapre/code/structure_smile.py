from rdkit import Chem
import csv

sdf_file = "structures.sdf"

# Initialize RDKit's SDF supplier to read the SDF file
supplier = Chem.SDMolSupplier(sdf_file, sanitize=False)
molecule_data = []  # This will store our database_id and SMILES pairs

for mol in supplier:
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)  # Try to sanitize the molecule
            # Extract DATABASE_ID, assuming it exists. If not, use a placeholder like 'No ID'
            database_id = mol.GetProp("DATABASE_ID") if mol.HasProp("DATABASE_ID") else "No ID"
            # Convert the molecule to a SMILES string
            smiles = Chem.MolToSmiles(mol)
            # Append the database_id and SMILES to our list
            molecule_data.append([database_id, smiles])
            print("Number of molecules loaded:" , len(molecule_data))
        except:
            print("Skipping a molecule due to sanitization error.")

# Print the number of molecules processed
print("Number of molecules processed:", len(molecule_data))

# Now, save the data to a CSV file
csv_file = "molecules.csv"
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["DATABASE_ID", "SMILES"])  # Writing the header
    writer.writerows(molecule_data)

print(f"Data saved to {csv_file}")