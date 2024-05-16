import matplotlib.pyplot as plt
import pandas as pd

file_path = 'Training_info/training_info.csv'
data = pd.read_csv(file_path)

plt.figure(figsize=(10, 8))  # Set the figure size

# Plot each metric
#plt.plot(data['Epoch'], data['Loglike'], label='Log-Likelihood', marker='o')
#plt.plot(data['Epoch'], data['NLL_x'], label='NLL X', marker='o')
#plt.plot(data['Epoch'], data['NLL_adj'], label='NLL Adj', marker='o')
plt.plot(data['Epoch'], data['KL_loss'], label='KL Divergence', marker='o')

# Adding labels and title
plt.xlabel('Epochs')
plt.ylabel('Metric Values')
plt.title('Training Metrics Over Epochs')

# Adding a legend
plt.legend()

# Show grid
plt.grid(True)

# Display the plot
plt.show()