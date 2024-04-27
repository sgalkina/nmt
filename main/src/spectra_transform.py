import os
import pandas as pd
import numpy as np
import dataset as dt

spec_max_mz = 2500
max_num_peaks = 100
min_intensity = 0.1

def process_spectrum(file_path):
    """ Process a single spectrum file and return the transformed data. """
    try:
        # spectrum = pd.read_csv(file_path, sep='\t', header=None)
        spectrum = pd.read_csv(file_path, sep='\t', header=None, names=['m/z', 'intensity'])
        spectrum = dt.FilterPeaks(spec_max_mz, min_intensity)(spectrum)
        spectrum = dt.Normalize(intensity=True, mass=True, rescale_intensity=True, max_mz=spec_max_mz)(spectrum)
        spectrum = dt.TopNPeaks(100)(spectrum)
        spectrum = dt.ToMZIntConcatAlt(max_num_peaks)(spectrum)

        if np.all(spectrum == 0):
            print(f"All zeros in file: {file_path}")
            return np.array([])  
        
        return spectrum
    
    except Exception as e:
        print(f"Failed to process {file_path}: {str(e)}")
        return np.array([])  

def process_all_spectra(directory):
    """ Process all spectra in the given directory and save to a CSV file. """
    transformed_data = []
    number_of_files_transformed = 0
    number_of_files_failed = 0
    filenames = [f for f in os.listdir(directory) if f.endswith('.txt') or f.endswith('.tsv')]
    for filename in filenames:
        file_path = os.path.join(directory, filename)
        transformed_spectrum = process_spectrum(file_path)
        ## Adding results with elements
        if transformed_spectrum.size > 0: 
            transformed_data.append(transformed_spectrum)
            number_of_files_transformed += 1
            print("Files transformed: ", number_of_files_transformed)
        else: 
            number_of_files_failed += 1
            print("Files failed: ", number_of_files_failed)

    # Convert list of arrays to a single numpy array
    transformed_data = np.vstack(transformed_data)
    # Save the transformed data to CSV
    transformed_file_path = '../data/transformed_spectra_for_all_spectra.csv'
    np.savetxt(transformed_file_path, transformed_data, delimiter=',', fmt='%s')
    print(f"Transformed data saved to {transformed_file_path}")

# Specify the directory containing the spectra files
directory = '../data/hmdb_experimental_msms_peak_lists'
process_all_spectra(directory)

transformed_data_path = '../data/transformed_spectra_for_all_spectra.csv'
df = pd.read_csv(transformed_data_path, header=None) 
print("Shape of the transformed data:", df.shape)

# Select the first 20 rows
first_20 = df.iloc[:20] 

new_file_path = '../data/transformed20.csv'
first_20.to_csv(new_file_path, index=False, header=False) 

print(f"First 20 columns saved to {new_file_path}")