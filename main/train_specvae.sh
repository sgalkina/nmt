#!/bin/bash
#SBATCH -p gpu                          # The partition/queue to run on, 'gpu' is standard
#SBATCH --gres=gpu:1                    # Request GPU resources
#SBATCH --job-name=spectrum_encoder     # Name of the job
#SBATCH --ntasks=1                      # Number of tasks (usually should be 1 for non-parallel jobs)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=250G                      # Memory needed (adjust if necessary)
#SBATCH --time=10-00:00:00              # Expected runtime (max)
#SBATCH --output=result_%j.out          # Standard output and error log

# Load necessary modules
module load cuda/11.8
module load gcc/11.2.0
module load anaconda3/2023.03-py3.10

# Activate conda environment
source ~/.bashrc
conda activate moflow

# Navigate to the directory
cd Multi-SpecVAE/main/src

# Train the spectrum encoder
python train_example.py 
