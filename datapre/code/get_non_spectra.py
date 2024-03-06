from rdkit import Chem
import pandas as pd

# Load the SDF file
sdf_path = 'structures.sdf'
sdf_supplier = Chem.SDMolSupplier(sdf_path)

# Load the xml_smiles.csv file
xml_smiles_path = 'xml_smiles.csv'
xml_smiles_df = pd.read_csv(xml_smiles_path)

# Load the existing zinc250k.csv file
zinc250k_path = 'zinc250k.csv'
zinc250k_df = pd.read_csv(zinc250k_path)

# Extract SMILES and potentially other identifiers if needed
sdf_smiles = [Chem.MolToSmiles(mol) for mol in sdf_supplier if mol is not None]

# Create a DataFrame for the SDF data
sdf_df = pd.DataFrame(sdf_smiles, columns=['SMILES'])

# Find the SMILES in the SDF data that are not in the XML SMILES data
unique_sdf_smiles = sdf_df[~sdf_df['SMILES'].isin(xml_smiles_df['SMILES'])]

# Add the new unique SMILES to the zinc250k DataFrame
unique_sdf_smiles['logP'] = 0
unique_sdf_smiles['qed'] = 0
unique_sdf_smiles['SAS'] = 0

# Rename the column to match the zinc250k dataset
unique_sdf_smiles = unique_sdf_smiles.rename(columns={'SMILES': 'smiles'})

# Append the new entries to the zinc250k DataFrame
updated_zinc250k_df = pd.concat([zinc250k_df, unique_sdf_smiles], ignore_index=True)

# Save the updated DataFrame back to zinc250k.csv
updated_zinc250k_df.to_csv('zinc250k2024.csv', index=False)

print("Number of non-spectra molecules added: ", len(unique_sdf_smiles))



"""
SDF FIL FORMAT OG INDHOLD: 

'DATABASE_ID', 'DATABASE_NAME', 'SMILES', 'INCHI_IDENTIFIER',
       'INCHI_KEY', 'FORMULA', 'MOLECULAR_WEIGHT', 'EXACT_MASS',
       'JCHEM_ACCEPTOR_COUNT', 'JCHEM_ATOM_COUNT',
       'JCHEM_AVERAGE_POLARIZABILITY', 'JCHEM_BIOAVAILABILITY',
       'JCHEM_DONOR_COUNT', 'JCHEM_FORMAL_CHARGE', 'JCHEM_GHOSE_FILTER',
       'JCHEM_IUPAC', 'ALOGPS_LOGP', 'JCHEM_LOGP', 'ALOGPS_LOGS',
       'JCHEM_MDDR_LIKE_RULE', 'JCHEM_NUMBER_OF_RINGS',
       'JCHEM_PHYSIOLOGICAL_CHARGE', 'JCHEM_PKA_STRONGEST_ACIDIC',
       'JCHEM_PKA_STRONGEST_BASIC', 'JCHEM_POLAR_SURFACE_AREA',
       'JCHEM_REFRACTIVITY', 'JCHEM_ROTATABLE_BOND_COUNT',
       'JCHEM_RULE_OF_FIVE', 'ALOGPS_SOLUBILITY', 'JCHEM_TRADITIONAL_IUPAC',
       'JCHEM_VEBER_RULE', 'HMDB_ID', 'GENERIC_NAME', 'SYNONYMS', 'ID',
       'ROMol', 'JCHEM_PKA'],
"""