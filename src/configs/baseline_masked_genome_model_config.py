"""
File Name:          baseline_masked_genome_model_config.py
Project:            bioseq-learning

File Description:

"""
from math import sqrt
from typing import Optional, List, Dict, Any


# experiment name associated with this set of configurations
# this could be used for various bookkeeping purposes
_experiment_name: str = \
    'baseline masked genome model prediction with transformer'

# indicator experimental with (much) smaller training set
# and validation/test sets are the same as the training set
_dry_run: bool = False

# random seed and deterministic flag for reproducible results
_random_seed: int = 0
_deterministic_cudnn_flag: bool = True

# "preferred" GPUs for computation device specification
# to use CPU only, set to None or empty list []; otherwise, set to a list
# of integers representing preferred GPUs for this experiment
_preferred_gpu_list: Optional[List[int]] = [0, 1, 2, 3]
# flag for using multiple GPUs (nn.DataParallel) for this experiment
_multi_gpu_flag: bool = False

# Nvidia apex mixed-precision training
_nvidia_amp_opt: bool = True
_nvidia_amp_opt_level: str = 'O3'


# dataset and dataloader parameters
_vld_ratio: float = 0.1
_tst_ratio: float = 0.1
_seq_len: int = 2000
_num_masks: float = 0.15
_max_num_paddings: int = 500
_dataloader_batch_size: int = 32
_dataloader_num_workers: int = 20
_max_num_trn_batches_per_epoch: int = 10000
_max_num_vld_batches_per_epoch: int = 10000


# transformer and network module configurations
# the embedding dimension for each "word"(A, T, G, C, and <padding>)
# must be dividable by the number of attention heads
_emb_dim: int = 4
# boolean indicator for positional encoding
_pos_enc: bool = True
_pos_enc_dropout: float = 0.0
_pos_enc_emb_scale: float = sqrt(_emb_dim)
# the number of attention heads must be a factor of the embedding dimension
_xfmr_enc_layer_num_attn_heads: int = 2
_xfmr_enc_layer_feedforward_dim: int = 1024
_xfmr_enc_layer_activation: str = 'relu'
_xfmr_enc_layer_dropout: float = 0.0
_xfmr_enc_num_layers: int = 3
_xfmr_enc_norm: bool = True


# training configurations
_max_num_epochs: int = 500
_early_stopping_patience: int = 32
_optimizer: str = 'SGD'
_optimizer_kwargs: Dict[str, Any] = {
    'lr': 1e-4,
    'momentum': 0.9,
}
_lr_scheduler: str = 'CosineAnnealingWarmRestarts'
_lr_scheduler_kwargs: Dict[str, Any] = {
    'T_0': 16,
    'eta_min': 1e-6,
}
# logging configurations
_num_trn_logs: int = 10


# dictionary that maps names of each configuration to their object
# e.g. 'experiment_name': _experiment_name, 'random_state': _random_state, etc.
# note that the local configuration variable names must start with '_', but
# the underscores are stripped away in the CONFIG dictionary
config: Dict[str, object] = {
    variable_name[1:]: variable
    for variable_name, variable in locals().items() if
    variable_name.startswith('_') and not variable_name.startswith('__')
}
