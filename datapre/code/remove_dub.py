import pandas as pd

# Load your dataset (replace 'your_dataset.csv' with the path to your dataset file)
df = pd.read_csv('zinc250k&non_spectra.csv')

# Remove duplicate rows based on the 'smiles' column
df_unique = df.drop_duplicates(subset='smiles', keep='first')

# Save the result to a new CSV file
df_unique.to_csv('zinc250k&non_spectra_unique.csv', index=False)

print(f"Original dataset size: {len(df)}, Unique dataset size: {len(df_unique)}")