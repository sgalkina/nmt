import os
from xml.etree import ElementTree as ET
import csv

#####################################################################
### EXTRACTS THE DATABSE_ID AND STRUCTURE_ID FOR ALL XML FILES

import os
from xml.etree import ElementTree as ET
import csv

def extract_data_from_xml(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Extract database_id
    database_id = root.find('.//database-id')
    database_id_text = database_id.text if database_id is not None else None

    # Check for MS-MS spectra
    ms_ms_peaks = root.findall('.//ms-ms-peak')
    has_spectra = len(ms_ms_peaks) > 0 and all(peak.get('nil') != 'true' for peak in ms_ms_peaks)

    return database_id_text, has_spectra

xml_folder_path = '/Users/nikolaikjaernielsen/Desktop/hmdb_experimental_msms_spectra/'

data = []

for filename in os.listdir(xml_folder_path):
    if filename.endswith('.xml'):
        full_path = os.path.join(xml_folder_path, filename)
        # Extract data
        database_id, has_spectra = extract_data_from_xml(full_path)
        data.append([database_id, has_spectra])

csv_file_path = '/Users/nikolaikjaernielsen/Desktop/database_id_spectra.csv'

# Write data to CSV
with open(csv_file_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['DATABASE_ID', 'HAS_SPECTRA'])
    csvwriter.writerows(data)

print(f"Data saved to {csv_file_path}")