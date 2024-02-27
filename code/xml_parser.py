import os
from xml.etree import ElementTree as ET
import csv

#####################################################################
### EXTRACTS THE DATABSE_ID AND STRUCTURE_ID FOR ALL XML FILES

def extract_ids_from_xml(xml_path):
    # Extract structure and database IDs from an XML file.
    tree = ET.parse(xml_path)
    root = tree.getroot()
    structure_id = root.find('.//structure-id')
    database_id = root.find('.//database-id')
    return structure_id.text if structure_id is not None else None, database_id.text if database_id is not None else None

# Path to your folder containing XML files
xml_folder_path = '/Users/nikolaikjaernielsen/Desktop/hmdb_experimental_msms_spectra/'

all_structure_ids = []
all_database_ids = []

# Loop through all files in the specified directory
for filename in os.listdir(xml_folder_path):
    if filename.endswith('.xml'):
        full_path = os.path.join(xml_folder_path, filename)
        # Extract IDs
        structure_id, database_id = extract_ids_from_xml(full_path)
        if structure_id and database_id:  # Ensure both IDs are not None
            all_structure_ids.append(structure_id)
            all_database_ids.append(database_id)

csv_file_path = '/Users/nikolaikjaernielsen/Desktop/xml_database_ids.csv'

# Open the CSV file in write mode
with open(csv_file_path, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csvwriter = csv.writer(csvfile)
    
    # Write the header
    csvwriter.writerow(['DATABASE_ID'])
    
    # Write each database ID as a new row
    for database_id in all_database_ids:
        csvwriter.writerow([database_id])

print(f"Database IDs saved to {csv_file_path}")
