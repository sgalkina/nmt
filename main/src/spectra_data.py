import os
import pandas as pd

directory = '/Users/nikolaikjaernielsen/Desktop/hmdb_experimental_msms_peak_lists/'
output_csv_path = '/Users/nikolaikjaernielsen/Desktop/spectra.csv'
spectra_list = []

# Define limit for the number of spectra to process (-1 for no limit)
limit = -1

file_paths = [os.path.join(directory, f) for f in os.listdir(directory)]
files_loaded = 0

for file_path in file_paths:
    with open(file_path, 'r') as file:
        spectrum = []  
        for line in file:
            try:
                mz, intensity = map(float, line.split())
                spectrum.append([mz, intensity])
            except ValueError:
            # Skip lines that can't be converted to floats
                continue
        files_loaded += 1
        print(files_loaded)
        spectra_list.extend(spectrum)
        # Break the loop if the limit is reached
        if limit > 0 and len(spectra_list) >= limit:
            break

spectra_df = pd.DataFrame(spectra_list, columns=['m/z', 'intensity'])

# Save the DataFrame to a CSV file
spectra_df.to_csv(output_csv_path, index=False)

print(f"Spectra data saved to {output_csv_path}")
