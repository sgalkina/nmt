import numpy as np
from collections import OrderedDict
from model import BaseModel

import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision as tv

class BaseVAE(BaseModel):
    def __init__(self, config, device=None):
        super().__init__(config, device)
        self.name = self.get_attribute('name', required=False, default='vae')

    def build_layers(self):
        self.encoder = None
        self.fc_mean = None
        self.fc_logvar = None
        self.decoder = None
        self.loss = None
        raise NotImplementedError("Build model method 'build_layers' is not implemented")

    def reparameterize(self, latent_dist):
        mean, logvar = latent_dist
        if self.training:
            z = self.sample(latent_dist)
        else:
            z = mean
        return z

    def encode(self, x):
        x = self.encoder(x)
        mu = self.fc_mean(x)
        logvar = self.fc_logvar(x)
        return mu, logvar

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        x, _, __ = self.forward_(x)
        return x

    def forward_(self, x):
        latent_dist = self.encode(x)
        z = self.reparameterize(latent_dist)
        x = self.decode(z)
        return x, z, latent_dist

    def get_layer_string(self):
        layers = self.layer_config
        layer_string = '-'.join(str(x) for x in layers[0]) + '-' + '-'.join(str(x) for x in np.array(layers[1])[1:])
        return layer_string

class SpecVEA(BaseVAE):
    def __init__(self, config, device=None):
        super(SpecVEA, self).__init__(config, device)
        self.name = self.get_attribute('name', required=False, default='specvae')
        self.beta = self.get_attribute('beta', required=False, default=1.0)
        self.build_layers()
        if self.device:
            self.to(self.device)

    def build_layers(self):
        # Build model layers
        self.layer_config = self.config['layer_config']
        self.encoder_layer_config = self.layer_config[0]
        self.decoder_layer_config = self.layer_config[1]

        # Encoder layers:
        encoder_layers_num = len(self.encoder_layer_config)
        encoder_layers = []
        for i in range(1, encoder_layers_num - 1):
            in_dim, out_dim = self.encoder_layer_config[i - 1], self.encoder_layer_config[i]
            encoder_layers.append(('en_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            encoder_layers.append(('en_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            encoder_layers.append(('en_act_%d' % i, nn.ReLU()))
        self.encoder = nn.Sequential(OrderedDict(encoder_layers))

        # Latent space layer (mu & log_var):
        in_dim, out_dim = \
            self.encoder_layer_config[encoder_layers_num - 2], \
            self.encoder_layer_config[encoder_layers_num - 1]
        self.latent_dim = out_dim
        self.fc_mean = nn.Linear(in_dim, out_dim)
        self.mean_batchnorm = nn.BatchNorm1d(out_dim)
        self.fc_logvar = nn.Linear(in_dim, out_dim)
        self.logvar_batchnorm = nn.BatchNorm1d(out_dim)

        # Sample from N(0., 1.)
        self.sample = SampleZ()

        # Decoder layers:
        decoder_layers_num = len(self.decoder_layer_config)
        decoder_layers = []
        for i in range(1, decoder_layers_num - 1):
            in_dim, out_dim = self.decoder_layer_config[i - 1], self.decoder_layer_config[i]
            decoder_layers.append(('de_lin_%d' % i, nn.Linear(in_dim, out_dim)))
            decoder_layers.append(('de_lin_batchnorm_%d' % i, nn.BatchNorm1d(out_dim)))
            decoder_layers.append(('de_act_%d' % i, nn.ReLU()))

        # Last layer of decoder:
        in_dim, out_dim = \
            self.decoder_layer_config[decoder_layers_num - 2], \
            self.decoder_layer_config[decoder_layers_num - 1]
        decoder_layers.append(('de_lin_%d' % (decoder_layers_num - 1), nn.Linear(in_dim, out_dim)))
        decoder_layers.append(('de_act_%d' % (decoder_layers_num - 1), nn.Sigmoid()))
        self.decoder = nn.Sequential(OrderedDict(decoder_layers))

        # Loss:
        self.input_size = self.encoder_layer_config[0]
        self.loss = VAELoss(self.beta)

    def encode(self, x):
        x = self.encoder(x)
        mu = self.mean_batchnorm(self.fc_mean(x))
        log_var = self.logvar_batchnorm(self.fc_logvar(x))
        return mu, log_var

class SampleZ(nn.Module):
    def forward(self, x):
        mu, log_sigma = x
        std = torch.exp(0.5 * log_sigma).to(mu.device)
        with torch.no_grad():
            epsilon = torch.randn_like(std).to(mu.device)
        return mu + std * epsilon

class VAELoss(nn.Module):
    def __init__(self, beta=1.0):
        super().__init__()
        self.beta = beta

    def forward(self, input, target, latent_dist):
        loss, _, __ = self.forward_(input, target, latent_dist)
        return loss

    def forward_(self, input, target, latent_dist):
        mean, logvar = latent_dist
        if torch.isnan(input).any():
            raise ValueError("input for VAELoss contains nan elements")
        # E[log P(X|z)]
        recon = torch.sum(F.binary_cross_entropy(input, target, reduction='none'), dim=1).mean()
        # D_KL(Q(z|X) || P(z|X))
        kld = self.beta * 0.5 * torch.sum(torch.exp(logvar) + torch.square(mean) - 1. - logvar, dim=1).mean()
        loss = recon + kld
        return loss, recon, kld