# GF-VAE
A Flow-based Variational Autoencoder for Molecule Generation

the code can be separated into two parts, skeleton of the flow part is from paper MoFlow, which you may approach at https://github.com/calvin-zcx/moflow, thanks for his great job. The detailed code is shown in model folder.

I used this to train: 

python train_model.py  --data_name hmdb  --batch_size  256  -learning_rate 0.0000 1--max_epochs 100 --gpu 0  --debug True  --save_dir=results/hmdb_2024   --b_n_flow 10  --b_hidden_ch 512,512  --a_n_flow 27  --a_hidden_gnn 64  --a_hidden_lin  512,64   --mask_row_size_list 1 --mask_row_stride_list 1  --noise_scale 0.6  --b_conv_lu 1

python train_model.py --data_name hmdb  --batch_size 128  --learning_rate 0.0001 --max_epochs 200 --gpu 0  --debug True  --save_dir=results/hmdb_2024_2  --b_n_flow 10  --b_hidden_ch 128  --a_n_flow 3 --a_hidden_gnn 64  --a_hidden_lin 64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 1 

python train_model.py  --data_name hmdb  --batch_size  256  --max_epochs 100 --gpu 0  --debug True  --save_dir=results/hmdb_2024   --b_n_flow 10  --b_hidden_ch 512,512  --a_n_flow 27  --a_hidden_gnn 64  --a_hidden_lin  512,64   --mask_row_size_list 1 --mask_row_stride_list 1  --noise_scale 0.6  --b_conv_lu 1

python train_model.py --data_name hmdb  --batch_size 32  --learning_rate 0.0001 --max_epochs 200 --gpu 0  --debug True  --save_dir=results/hmdb_2024_2  --b_n_flow 10  --b_hidden_ch 128,128  --a_n_flow 27 --a_hidden_gnn 64  --a_hidden_lin 128,64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.3 --b_conv_lu 1 

python train_model.py --data_name hmdb  --batch_size 256  --learning_rate 0.001 --max_epochs 200 --gpu 0  --debug True  --save_dir=results/hmdb_2024_2  --b_n_flow 4  --b_hidden_ch 64,64  --a_n_flow 6 --a_hidden_gnn 64  --a_hidden_lin 64,32  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.3 --b_conv_lu 1 

python train_model.py --data_name hmdb  --batch_size 128  --learning_rate 0.00001 --max_epochs 200 --gpu 0  --debug True  --save_dir=results/hmdb_15-04 --b_n_flow 10  --b_hidden_ch 128  --a_n_flow 3 --a_hidden_gnn 64  --a_hidden_lin 64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 1 

And I use this to reconstruct: 

python generate.py --model_dir results/hmdb_2024_2 -snapshot model_snapshot_epoch_200 --gpu 0 --data_name hmdb --hyperparams-path moflow-params.json --batch-size 256 --reconstruct  2>&1 | tee hmdb_2024_april_reconstruct_results.txt


