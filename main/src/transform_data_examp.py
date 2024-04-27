import dataset as dt
import numpy as np
import pandas as pd
# from disentanglement import dataset as dt

spec_max_mz = 2500
max_num_peaks = 1000
min_intensity = 0.1

# 42.794 rows and 200 columns and keep track of filenames

file_path = '../data/hmdb_experimental_msms_peak_lists/HMDB0000001_msms_1_2_experimental.txt'
# file_path = '../data/spectra.csv'
spectrum = pd.read_csv(file_path, sep='\t')

# spectra_example = spectra_data[:,500]

spectrum = dt.FilterPeaks(max_mz=spec_max_mz, min_intensity=min_intensity)(spectrum)
print("FilterPeaks COMPLETED")

spectrum = dt.Normalize(intensity=True, mass=True, rescale_intensity=True, max_mz=spec_max_mz)(spectrum)
print("Normalize COMPLETED")
print("SPECTRUM AFTER NORMALIZE:", spectrum)

spectrum = dt.TopNPeaks(100)(spectrum)
print("TopNPeaks COMPLETED")
print("SPECTRUM AFTER TOPNPEAKS:", spectrum)

spectrum = dt.ToMZIntConcatAlt(max_num_peaks=max_num_peaks)(spectrum)
print("ToMZIntConcatAlt COMPLETED")




transformed_file_path = '../data/transformed_spectra_for_one_spectra.csv' 
np.savetxt(transformed_file_path, spectrum, delimiter=',', fmt='%s')
print(f"Transformed data saved to {transformed_file_path}")