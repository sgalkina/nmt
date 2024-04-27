import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

original_spectrum = pd.read_csv('../data/transformed_spectra.csv')
reconstructed_spectrum = pd.read_csv('../results/spectra_reconstructions.csv')

"""
first_20_originals = original_spectrum[:20]
first_20_reconstructed = reconstructed_spectrum[:20]

first_20_originals.to_csv('first_20_originals.csv', index=False)
first_20_reconstructed.to_csv('first_20_reconstructions.csv', index=False)""" 

print(original_spectrum.shape)
print(reconstructed_spectrum.shape)

num_plots = 2
start_value = 11
fig, axes = plt.subplots(num_plots, 1, figsize=(10, 15))
for i in range(num_plots):
    if isinstance(original_spectrum, pd.DataFrame):
        original_data = original_spectrum.iloc[i+start_value].values 
    else:
        original_data = original_spectrum[i]  

    if isinstance(reconstructed_spectrum, pd.DataFrame):
        reconstructed_data = reconstructed_spectrum.iloc[i+start_value].values  
    else:
        reconstructed_data = reconstructed_spectrum[i]

    axes[i].plot(original_data, label='Original')
    axes[i].plot(reconstructed_data, label='Reconstructed', linestyle='--')
    axes[i].set_title(f'Spectrum {i + start_value +1}')
    axes[i].legend()

plt.tight_layout()
plt.show()