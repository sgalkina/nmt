import json
import os
import sys
# for linux env.
sys.path.insert(0,'..')
import argparse
from distutils.util import strtobool
import torch
import torch.nn as nn
import numpy as np

from data.data_loader import NumpyTupleDataset
from mflow.models.hyperparams import Hyperparameters
# from mflow.models.model import MoFlow, rescale_adj
from mflow.models.graphvae_flow import MoFlow, rescale_adj
from mflow.models.utils import check_validity, save_mol_png
from torch.autograd import Variable 

import time
from mflow.utils.timereport import TimeReport
from mflow.generate import generate_mols
import torch.nn.functional as F

from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')     

import functools
print = functools.partial(print, flush=True)


def get_parser():
    parser = argparse.ArgumentParser()
    # data I/O
    parser.add_argument('-i', '--data_dir', type=str, default='data', help='Location for the dataset')
    parser.add_argument('--data_name', type=str, default='qm9', choices=['qm9', 'zinc250k'], help='dataset name')
    # parser.add_argument('-f', '--data_file', type=str, default='qm9_relgcn_kekulized_ggnp.npz', help='Name of the dataset')
    parser.add_argument('-o', '--save_dir', type=str, default='results/qm9',
                        help='Location for parameter checkpoints and samples')
    parser.add_argument('-t', '--save_interval', type=int, default=20,
                        help='Every how many epochs to write checkpoint/samples?')
    parser.add_argument('-r', '--load_params', type=int, default=0,
                        help='Restore training from previous model checkpoint? 1 = Yes, 0 = No')
    parser.add_argument('--load_snapshot', type=str, default='', help='load the model from this path')
    # optimization
    parser.add_argument('-l', '--learning_rate', type=float, default=0.001, help='Base learning rate')
    parser.add_argument('-e', '--lr_decay', type=float, default=0.999995,
                        help='Learning rate decay, applied every step of the optimization')
    parser.add_argument('-x', '--max_epochs', type=int, default=5000, help='How many epochs to run in total?')
    parser.add_argument('-g', '--gpu', type=int, default=0, help='GPU Id to use')
    parser.add_argument('--save_epochs', type=int, default=1, help='in how many epochs, a snapshot of the model'
                                                                   ' needs to be saved?')
    # data loader
    parser.add_argument('-b', '--batch_size', type=int, default=12, help='Batch size during training per GPU')
    parser.add_argument('--shuffle', type=strtobool, default='false', help='Shuffle the data batch')
    parser.add_argument('--num_workers', type=int, default=0, help='Number of workers in the data loader')

    # # evaluation
    # parser.add_argument('--sample_batch_size', type=int, default=16,
    #                     help='How many samples to process in paralell during sampling?')
    # reproducibility
    
    # For encoder & decoder model
    parser.add_argument('--enc_conv_dim', default=[[1024,512,256,128],128], help='number of conv filters in the grpah encoder')
    parser.add_argument('--enc_linear_dim', default=[1024,512,256,128], help='linear dims in graph encoder')
    parser.add_argument('--dec_dim', default=[128, 128], help='linear dims in graph decoder')
    # For bonds
    parser.add_argument('--b_n_flow', type=int, default=4,
                        help='Number of masked glow coupling layers per block for bond tensor')
    parser.add_argument('--b_n_block', type=int, default=1, help='Number of glow blocks for bond tensor')
    parser.add_argument('--b_hidden_ch', type=str, default="128,128",
                        help='Hidden channel list for bonds tensor, delimited list input ')
    parser.add_argument('--b_conv_lu', type=int, default=1, choices=[0, 1, 2],
                        help='0: InvConv2d for 1*1 conv, 1:InvConv2dLU for 1*1 conv, 2: No 1*1 conv, '
                             'swap updating in the coupling layer')
    # For atoms
    parser.add_argument('--a_n_flow', type=int, default=4,
                        help='Number of masked flow coupling layers per block for atom matrix')
    parser.add_argument('--a_n_block', type=int, default=1, help='Number of flow blocks for atom matrix')
    parser.add_argument('--a_hidden_gnn', type=str, default="64,",
                        help='Hidden dimension list for graph convolution for atoms matrix, delimited list input ')
    parser.add_argument('--a_hidden_lin', type=str, default="128,64",
                        help='Hidden dimension list for linear transformation for atoms, delimited list input ')
    parser.add_argument('--mask_row_size_list', type=str, default="1,",
                        help='Mask row size list for atom matrix, delimited list input ')
    parser.add_argument('--mask_row_stride_list', type=str, default="1,",
                        help='Mask row stride list for atom matrix, delimited list input')
    # General
    parser.add_argument('-s', '--seed', type=int, default=1, help='Random seed to use')
    parser.add_argument('--debug', type=strtobool, default='true', help='To run training with more information')
    parser.add_argument('--learn_dist', type=strtobool, default='true', help='learn the distribution of feature matrix')
    parser.add_argument('--noise_scale', type=float, default=0.6, help='x + torch.rand(x.shape) * noise_scale')

    return parser


def get_kl_loss(mu, logvar):
    kld_loss = torch.mean(-0.5 * torch.sum(1 + logvar - mu ** 2 - logvar.exp(), dim=1), dim=0)
    return kld_loss


def train():
    start = time.time()
    print("Start at Time: {}".format(time.ctime()))
    parser = get_parser()
    args = parser.parse_args()

    # Device configuration
    device = -1
    multigpu = False
    if args.gpu >= 0:
        # signle gpu
        # device = args.gpu
        device = torch.device('cuda:'+str(args.gpu) if torch.cuda.is_available() else 'cpu')
    elif args.gpu == -1:
        # cpu
        device = torch.device('cpu')
    else:
        # multigpu, can be slower than using just 1 gpu
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        multigpu = True

    debug = args.debug

    print('input args:\n', json.dumps(vars(args), indent=4, separators=(',', ':')))  # pretty print args

    # Model configuration
    enc_conv_dim = [[256, 128],128]
    enc_linear_dim = [128]
    dec_dim = [128, 128]

    # enc_conv_dim = [[1024,512,256,128],128]
    # enc_linear_dim = [1024,512,256,128]
    # dec_dim = [128, 128]

#     enc_conv_dim = [[16],16]
#     enc_linear_dim = [16]
#     dec_dim = [128]
    
    
    b_hidden_ch = [int(d) for d in args.b_hidden_ch.strip(',').split(',')]
    a_hidden_gnn = [int(d) for d in args.a_hidden_gnn.strip(',').split(',')]
    a_hidden_lin = [int(d) for d in args.a_hidden_lin.strip(',').split(',')]
    mask_row_size_list = [int(d) for d in args.mask_row_size_list.strip(',').split(',')]
    mask_row_stride_list = [int(d) for d in args.mask_row_stride_list.strip(',').split(',')]
    if args.data_name == 'qm9':
        from data import transform_qm9
        data_file = 'qm9_relgcn_kekulized_ggnp.npz'
        transform_fn = transform_qm9.transform_fn
        atomic_num_list = [6, 7, 8, 9, 0]
        # atomic_num_list = [6, 7, 8, 9, 15, 16, 17, 34, 35, 53, 0]
        b_n_type = 4
        b_n_squeeze = 3    # 3
        a_n_node = 9
        a_n_type = len(atomic_num_list)  # 5
        valid_idx = transform_qm9.get_val_ids()  # len: 13,082, total data: 133,885
    elif args.data_name == 'zinc250k':
        from data import transform_zinc250k
        data_file = 'zinc250k_relgcn_kekulized_ggnp.npz'
        transform_fn = transform_zinc250k.transform_fn_zinc250k
        atomic_num_list = transform_zinc250k.zinc250_atomic_num_list  # [6, 7, 8, 9, 15, 16, 17, 35, 53, 0]
        # mlp_channels = [1024, 512]
        # gnn_channels = {'gcn': [16, 128], 'hidden': [256, 64]}
        b_n_type = 4
        b_n_squeeze = 19   # 19
        a_n_node = 38
        a_n_type = len(atomic_num_list)  # 10
        valid_idx = transform_zinc250k.get_val_ids()
    else:
        raise ValueError('Only support qm9 and zinc250k right now. '
                         'Parameters need change a little bit for other dataset.')

    model_params = Hyperparameters(b_n_type=b_n_type,  # 4,
                                   b_n_flow=args.b_n_flow,
                                   b_n_block=args.b_n_block,
                                   b_n_squeeze=b_n_squeeze,
                                   b_hidden_ch=b_hidden_ch,
                                   b_affine=True,
                                   b_conv_lu=args.b_conv_lu,
                                   a_n_node=a_n_node,
                                   a_n_type=a_n_type,
                                   a_hidden_gnn=a_hidden_gnn,
                                   a_hidden_lin=a_hidden_lin,
                                   a_n_flow=args.a_n_flow,
                                   a_n_block=args.a_n_block,
                                   mask_row_size_list=mask_row_size_list,
                                   mask_row_stride_list=mask_row_stride_list,
                                   a_affine=True,
                                   device = args.gpu,
                                   learn_dist=args.learn_dist,
                                   seed=args.seed,
                                   noise_scale=args.noise_scale,
                                   enc_conv_dim = enc_conv_dim,
                                   enc_linear_dim = enc_linear_dim,
                                   dec_dim = dec_dim
                                   )
    print('Model params:')
    model_params.print()
    model = MoFlow(model_params)
    os.makedirs(args.save_dir, exist_ok=True)
    model.save_hyperparams(os.path.join(args.save_dir, 'moflow-params.json'))
    if torch.cuda.device_count() > 1 and multigpu:
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
        model = nn.DataParallel(model)
    else:
        multigpu = False
    model = model.to(device)

    # Datasets:
    dataset = NumpyTupleDataset.load(os.path.join(args.data_dir, data_file), transform=transform_fn)  # 133885

#     if len(valid_idx) > 0:
#         train_idx = [t for t in range(len(dataset)) if t not in valid_idx]  # 120803 = 133885-13082
#         # n_train = len(train_idx)  # 120803
#         train = torch.utils.data.Subset(dataset, train_idx)  # 120,803
#         test = torch.utils.data.Subset(dataset, valid_idx)  # 13,082
#     else:
#         torch.manual_seed(args.seed)
#         train, test = torch.utils.data.random_split(
#             dataset,
#             [int(len(dataset) * 0.8), len(dataset) - int(len(dataset) * 0.8)])

    train_dataloader = torch.utils.data.DataLoader(dataset, batch_size=args.batch_size,
                                                   shuffle=args.shuffle, num_workers=args.num_workers)

    print('==========================================')
    print('Load data done! Time {:.2f} seconds'.format(time.time() - start))
    print('Data shuffle: {}, Number of data loader workers: {}!'.format(args.shuffle, args.num_workers))
    if args.gpu >= 0:
        print('Using GPU device:{}!'.format(args.gpu))
#     print('Num Train-size: {}'.format(len(train)))
    print('Num Minibatch-size: {}'.format(args.batch_size))
    print('Num Iter/Epoch: {}'.format(len(train_dataloader)))
    print('Num epoch: {}'.format(args.max_epochs))
    print('==========================================')

    # Loss and optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)

    # Train the models
    iter_per_epoch = len(train_dataloader)
    log_step = args.save_interval  # 20 default
    tr = TimeReport(total_iter=args.max_epochs * iter_per_epoch)
    uniques = []
    for epoch in range(args.max_epochs):
        print("In epoch {}, Time: {}".format(epoch+1, time.ctime()))
        for i, batch in enumerate(train_dataloader):
            optimizer.zero_grad()
            # turn off shuffle to see the order with original code
            x = batch[0].to(device)  # (256,9,5)
            adj = batch[1].to(device)  # (256,4,9, 9)
            adj_normalized = rescale_adj(adj).to(device)

            # TODO: the part of the code that replaces the spectrum
            adj_flatten = torch.flatten(adj, start_dim=1)
            x_flatten = torch.flatten(x, start_dim=1)
            context = torch.cat([adj_flatten, x_flatten], 1)
            _, L = context.shape

            context_flatten = Variable(F.pad(context, (1, 9*5*10-1-L), "constant", 0), requires_grad=True)

            # Forward, backward and optimize
            z, sum_log_det_jacs, vae_para = model(adj, x, adj_normalized, context=context_flatten)
            kl_loss = get_kl_loss(vae_para[0], vae_para[1])
            z_plain = torch.cat((z[0].reshape(z[0].shape[0],-1),z[1].reshape(z[1].shape[0],-1)),dim=1)
            recon_loss = torch.sqrt(nn.MSELoss()(z_plain, vae_para[2]))

            adjm, xm = model.reverse(vae_para[2], context=context_flatten)
            recon_loss_adj = torch.sqrt(nn.MSELoss()(adj, adjm))
            recon_loss_x = torch.sqrt(nn.MSELoss()(x, xm))

            if multigpu:
                nll = model.module.log_prob(z, sum_log_det_jacs, vae_para)
            else:
                nll = model.log_prob(z, sum_log_det_jacs, vae_para)
            coef = 1000
            loss = nll[0]*coef + nll[1]*coef + kl_loss + recon_loss*coef + recon_loss_adj + recon_loss_x
            # loss = nll[0] + nll[1] + kl_loss + recon_loss
            loss.backward()
            optimizer.step()
            tr.update()

            # Print log info
            if (i+1) % log_step == 0:  # i % args.log_step == 0:
                print('Epoch [{}/{}], Iter [{}/{}], loglik: {:.5f}, nll_x: {:.5f},'
                      ' nll_adj: {:.5f}, kl_loss: {:.5f}, {:.2f} sec/iter, {:.2f} iters/sec:  '
                      ' recon_loss: {:.5f}, recon_loss_adj: {:.5f}, recon_loss_x: {:.5f}'.
                      format(epoch+1, args.max_epochs, i+1, iter_per_epoch,
                             loss.item(), nll[0].item(), nll[1].item(), kl_loss,
                             tr.get_avg_time_per_iter(), tr.get_avg_iter_per_sec(),
                             recon_loss, recon_loss_adj, recon_loss_x))
                tr.print_summary()

#         if debug:
#             def print_validity(ith):
#                 model.eval()
#                 if multigpu:
#                     adj, x = generate_mols(model.module, batch_size=100, device=device)
#                 else:
#                     adj, x = generate_mols(model, batch_size=1000, device=device)
#                 result = check_validity(adj, x, atomic_num_list, correct_validity=True)
#                 valid_mols = result['valid_mols']
#                 uniques.append(result['unique_ratio']*result['valid_ratio']/100.0)
# #                 mol_dir = os.path.join(args.save_dir, 'generated_{}'.format(ith))
# #                 os.makedirs(mol_dir, exist_ok=True)
# #                 for ind, mol in enumerate(valid_mols):
# #                     save_mol_png(mol, os.path.join(mol_dir, '{}.png'.format(ind)))
#                 model.train()
#             print_validity(epoch+1)

        # The same report for each epoch
        print('Epoch [{}/{}], Iter [{}/{}], loglik: {:.5f}, nll_x: {:.5f},'
              ' nll_adj: {:.5f}, kl_loss: {:.5f}, {:.2f} sec/iter, {:.2f} iters/sec '
                ' recon_loss: {:.5f}, recon_loss_adj: {:.5f}, recon_loss_x: {:.5f}'.
              format(epoch + 1, args.max_epochs, -1, iter_per_epoch,
                     loss.item(), nll[0].item(), nll[1].item(), kl_loss,
                     tr.get_avg_time_per_iter(), tr.get_avg_iter_per_sec(),
                     recon_loss.item(), recon_loss_adj.item(), recon_loss_x.item()))
        tr.print_summary()

        # Save the model checkpoints
        save_epochs = args.save_epochs
        if save_epochs == -1:
            save_epochs = args.max_epochs
        if (epoch + 1) % save_epochs == 0:
            if multigpu:
                torch.save(model.module.state_dict(), os.path.join(
                args.save_dir, 'model_snapshot_epoch_{}'.format(epoch + 1)))
            else:
                torch.save(model.state_dict(), os.path.join(
                args.save_dir, 'model_snapshot_epoch_{}'.format(epoch + 1)))
            tr.end()

#     print("[Training Ends], Start at {}, End at {}".format(time.ctime(start), time.ctime()))
#     print('totally {} seconds '.format(time.time()-start))
#     print('{} seconds per epoch'.format((time.time()-start)/args.max_epochs))
#     print(uniques)
#     with open(os.path.join(args.save_dir, 'unique.txt'),'w') as f:
#         for line in uniques:
#             f.write(str(round(line,3)))
#             f.write('\n')

#     model.load_state_dict(torch.load(os.path.join(args.save_dir, 'model_snapshot_epoch_50')))
#     model = model.to(device)
#     if args.data_name == 'qm9':
#         aggre_embeds = np.ones((1,9,128))
#     else:
#         aggre_embeds = np.ones((1,38,10))
#     for i, batch in enumerate(train_dataloader):
#         x = batch[0].to(device)  # (256,9,5)
#         adj = batch[1].to(device) # (256,4,9, 9)
#         atom_embed = model.graphenc.embed(adj,x)
#         aggre_embeds = np.concatenate((aggre_embeds, atom_embed.cpu().numpy()),axis=0)          
#     if adj.shape[2]==9:    
#         np.save('qm9node.npy',aggre_embeds[1:])
          


if __name__ == '__main__':
    # with torch.autograd.set_detect_anomaly(True):
    train()
# 42.4s/epoch for qm9_datatest_5_13
# 147.4s/epoch for zinc250k_datatest_5_15


# python train_model.py --data_name qm9  --batch_size 256  --max_epochs 50 --gpu 0  --debug True  --save_dir=results/qm9  --b_n_flow 5  --b_hidden_ch 128  --a_n_flow 13 --a_hidden_gnn 64  --a_hidden_lin 64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 1

# python train_model.py  --data_name zinc250k  --batch_size  256  --max_epochs 100 --gpu 0  --debug True  --save_dir=results/zinc250k_model_small_size/zinc_3_9  --b_n_flow 3  --b_hidden_ch 16  --a_n_flow 9 --a_hidden_gnn 16  --a_hidden_lin 16 --mask_row_size_list 1 --mask_row_stride_list 1  --noise_scale 0.6  --b_conv_lu 2

# python train_model.py  --data_name zinc250k  --batch_size  256  --max_epochs 100 --gpu 0  --debug True  --save_dir=results/zinc250k_datatest_5_15   --b_n_flow 5  --b_hidden_ch 256,128  --a_n_flow 15  --a_hidden_gnn 256  --a_hidden_lin  256,64   --mask_row_size_list 1 --mask_row_stride_list 1  --noise_scale 0.6  --b_conv_lu 2  

# python train_model.py --data_name qm9  --batch_size 256  --max_epochs 100 --gpu 0  --debug True  --save_dir=results/qm9_datatest_5_13  --b_n_flow 5  --b_hidden_ch 128  --a_n_flow 13 --a_hidden_gnn 64  --a_hidden_lin 64  --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 1 

# python train_model.py  --data_name zinc250k  --batch_size  256  --max_epochs 100 --gpu 1  --debug False  --save_dir=results/zinc250k_datatest_5_15   --b_n_flow 5  --b_hidden_ch 256,128  --a_n_flow 15  --a_hidden_gnn 256  --a_hidden_lin  256,64   --mask_row_size_list 1 --mask_row_stride_list 1 --noise_scale 0.6 --b_conv_lu 2


# python train_model.py  --data_name zinc250k  --batch_size  128  --max_epochs 100 --gpu 0  --debug True  --save_dir=results/zinc250k_test001_2squeeze   --b_n_flow 5  --b_hidden_ch 128,128  --a_n_flow 15  --a_hidden_gnn 256  --a_hidden_lin  512,64   --mask_row_size_list 1 --mask_row_stride_list 1 --b_conv_lu 2