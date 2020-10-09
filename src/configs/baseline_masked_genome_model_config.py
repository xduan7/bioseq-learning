"""
File Name:          baseline_masked_genome_model_config.py
Project:            bioseq-learning

File Description:

"""
import logging
from math import sqrt
from types import MappingProxyType
from typing import Optional, List, Dict, Union, Any


__LOGGER = logging.getLogger(__name__)


# experiment name associated with this set of configurations
# this could be used for various bookkeeping purposes
_experiment_name: str = \
    'baseline masked genome model prediction with transformer'

# indicator for Microsoft NNI hyper-parameter search
# note that the searching parameters from NNI would replace some of the
# parameters listed below in this file
_nni_search: bool = True

# indicator experimental with (much) smaller training set
# and validation/test sets are the same as the training set
_dry_run: bool = True

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
_nvidia_amp_opt: bool = False
_nvidia_amp_opt_level: str = 'O3'


# dataset and dataloader parameters
_vld_ratio: float = 0.01
_tst_ratio: float = 0.01
_seq_len: int = 2000
_num_masks: float = 0.10
_max_num_paddings: int = 0
_dataloader_batch_size: int = 32
_dataloader_num_workers: int = _dataloader_batch_size
_max_num_trn_batches_per_epoch: int = 10000
_max_num_vld_batches_per_epoch: int = 10000
_max_num_tst_batches: int = 10000


# transformer and network module configurations
# the embedding dimension for each "word"(A, T, G, C, and <padding>)
# must be dividable by the number of attention heads
_emb_dim: int = 1
# boolean indicator for positional encoding
_pos_enc: bool = True
_pos_enc_dropout: float = 0.0
_pos_enc_emb_scale: float = sqrt(_emb_dim)
# the number of attention heads must be a factor of the embedding dimension
_xfmr_enc_layer_num_attn_heads: int = 1
_xfmr_enc_layer_feedforward_dim: int = 1024
_xfmr_enc_layer_activation: str = 'relu'
_xfmr_enc_layer_dropout: float = 0.1
_xfmr_enc_num_layers: int = 1
_xfmr_enc_norm: bool = True


# training configurations
_max_num_epochs: int = 500
_early_stopping_patience: int = 32
# similar optimizer to the original BERT
# reference: https://arxiv.org/pdf/1706.03762.pdf
_optimizer: str = 'Adam'
_optimizer_kwargs: Dict[str, Any] = {
    'lr': 2e-3,
    'betas': (0.9, 0.98),
    'eps': 1e-9,
    # 'momentum': 0.9,
    # 'weight_decay': 1e-5,
}
_max_grad_norm: Union[int, float] = 10
_lr_scheduler: str = 'CosineAnnealingWarmRestarts'
_lr_scheduler_kwargs: Dict[str, Any] = {
    'T_0': 16,
    'eta_min': 1e-4,
}
# logging configurations
_num_trn_logs: int = 10

# adjust configurations if it's a dry run
if _dry_run:
    __dry_run_factor: int = 10
    _warning_msg = \
        f'reducing the data length (sequences and paddings), ' \
        f'the number of training and validation batches, and ' \
        f'the number of training and early stopping patience epochs ' \
        f'by a factor of {__dry_run_factor} for dry run ... '
    __LOGGER.warning(_warning_msg)

    _seq_len //= __dry_run_factor
    _max_num_paddings //= __dry_run_factor
    _max_num_trn_batches_per_epoch //= __dry_run_factor
    _max_num_vld_batches_per_epoch //= __dry_run_factor
    _max_num_epochs //= __dry_run_factor
    _early_stopping_patience //= __dry_run_factor


# dictionary that maps names of each configuration to their object
# e.g. 'experiment_name': _experiment_name, 'random_state': _random_state, etc.
# note that the local configuration variable names must start with '_', but
# the underscores are stripped away in the CONFIG dictionary
config: MappingProxyType = MappingProxyType({
    variable_name[1:]: variable
    for variable_name, variable in locals().items() if
    variable_name.startswith('_') and not variable_name.startswith('__')
})
