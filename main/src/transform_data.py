import dataset as dt
import numpy as np
import pandas as pd
# from disentanglement import dataset as dt
import torchvision.transforms as tv

spec_max_mz = 2500
max_num_peaks = 1000
min_intensity = 0.1
n_samples = -1 # -1 if all


train_data, valid_data, test_data = None, None, None
file_path = '../data/spectra.csv'
spectra_data = pd.read_csv(file_path)

# df_train, df_valid, df_test = MoNA.get_by_split(data_path, columns=['spectrum'])

# spectra_data = dt.SplitSpectrum()(spectra_data)
#print("SplitSpectrum COMPLETED")

spectra_data = dt.FilterPeaks(max_mz=spec_max_mz, min_intensity=min_intensity)(spectra_data)
print("FilterPeaks COMPLETED")

spectra_data = dt.Normalize(intensity=True, mass=True, max_mz=spec_max_mz)(spectra_data)
print("Normalize COMPLETED")

spectra_data = dt.ToMZIntConcatAlt(max_num_peaks=max_num_peaks)(spectra_data)
print("ToMZIntConcatAlt COMPLETED")

transformed_file_path = '../data/transformed_spectra.csv' 
np.savetxt(transformed_file_path, spectra_data, delimiter=',', fmt='%s')
print(f"Transformed data saved to {transformed_file_path}")