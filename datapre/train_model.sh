#!/bin/bash
#SBATCH -p gpu                      # The partition/queue to run on, 'gpu' is standard
#SBATCH --gres=gpu:1                # Request GPU resources
#SBATCH --job-name=moflow           # Name of the job
#SBATCH --ntasks=1                  # Number of tasks (usually should be 1 for non-parallel jobs)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=250G                  # Memory needed (adjust if necessary)
#SBATCH --time=6-00:00:00           # Expected runtime (max)
#SBATCH --output=result_%j.out      # Standard output and error log

# Load necessary modules
module load cuda/11.8
module load gcc/11.2.0
module load anaconda3/2023.03-py3.10

# Activate conda environment
source ~/.bashrc
conda activate moflow

# Navigate to the data directory and preprocess the data
cd moflow-master/data
python data_preprocess.py --data_name zinc250k

# Navigate to the mflow directory and run your machine learning training script
cd ../mflow
python train_model.py --data_name zinc250k  --batch_size 256  --max_epochs 50 --gpu 0  --debug True  --save_dir=results/qm9_64gnn_128-64lin_1-1mask_0d6noise_convlu1  --b_n_flow 10  --b_hidden_ch 128,128  --a_n_flow 27 --a_hidden_gnn 64  --a_hidden_lin 128,64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 1  2>&1 | tee qm9_64gnn_128-64lin_1-1mask_0d6noise_convlu1.log