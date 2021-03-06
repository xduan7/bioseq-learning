"""
File Name:          baseline_masked_genome_model_config.py
Project:            bioseq-learning

File Description:

"""
import logging
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
_nni_search: bool = False

# indicator experimental with (much) smaller training set
# and validation/test sets are the same as the training set
_dry_run: bool = False
_dry_run_contig_len: int = 10000
_dry_run_num_contigs: int = 10
_dry_run_num_ins_per_contig: int = 0
_dry_run_acc_target: float = 0.85

# random seed and deterministic flag for reproducible results
_random_seed: int = 0
_deterministic_cudnn_flag: bool = True

# "preferred" GPUs for computation device specification
# to use CPU only, set to None or empty list []; otherwise, set to a list
# of integers representing preferred GPUs for this experiment
_preferred_gpu_list: Optional[Union[List[int], str]] = \
    'all' if _nni_search else [0, 1, 2, 3, 4, 5, 6, 7]
# flag for using multiple GPUs (nn.DataParallel) for this experiment
_multi_gpu_flag: bool = True
# number of chunks for GPipe pipeline parallelization, only works when
# multi-GPU is enabled and necessary, must be a factor of batch size
_num_gpipe_chunks: int = 16

# Nvidia apex mixed-precision training
_nvidia_amp_opt: bool = False
_nvidia_amp_opt_level: str = 'O1'


# dataset and dataloader parameters
_trn_ratio: float = 1.00
_vld_ratio: float = 0.01
_tst_ratio: float = 0.01
_seq_len: int = 1000
_num_masks: float = 0.05
_max_num_paddings: Union[int, float] = 0.0
_dataloader_batch_size: int = 256
_dataloader_num_workers: int = _dataloader_batch_size
_max_num_trn_batches_per_epoch: int = 1000
_max_num_vld_batches_per_epoch: int = 1000
_max_num_tst_batches: int = 10000

# length of k-mer conversion for genomes
_kmer_len: int = 3

# transformer and network module configurations
# the embedding dimension for each "word"(A, T, G, C, and <padding>)
# - must be dividable by the number of attention heads
# - must be dividable by 2 with positional encoding
_emb_dim: int = 128
# boolean indicator for positional encoding
_pos_enc: bool = True
_pos_enc_dropout: float = 0.0
_pos_enc_trainable: bool = False
# _pos_enc_emb_scale: float = sqrt(_emb_dim)
# the number of attention heads must be a factor of the embedding dimension
_xfmr_enc_layer_num_attn_heads: int = max(_emb_dim // 32, 2)
_xfmr_enc_layer_feedforward_dim: int = 256
_xfmr_enc_layer_activation: str = 'relu'
# TODO: need to investigate xfmr dropout (no effect)
_xfmr_enc_layer_dropout: float = 0.0
_xfmr_enc_layer_norm: bool = False
_xfmr_enc_num_layers: int = 16
_xfmr_attn_mask: bool = False
_xfmr_padding_mask: bool = False
_xfmr_init_gain: float = 1.0 / (_xfmr_enc_num_layers ** 0.5)

# training configurations
_max_num_epochs: int = 10000
_early_stopping_patience: int = 500
# similar optimizer to the original BERT
# reference: https://arxiv.org/pdf/1706.03762.pdf
_optimizer: str = 'AdamW'
_optimizer_weight_decay: float = 1e-4
_optimizer_kwargs: Dict[str, Any] = {
    'lr': 1e-4,
    # 'amsgrad': True,
    # 'betas': (0.9, 0.98),
    # 'eps': 1e-9,
    # 'momentum': 0.9,
    'weight_decay': _optimizer_weight_decay,
}
_max_grad_norm: Optional[Union[int, float]] = 10.0
_warmup_factor: float = 0.01
_num_warmup_epochs: int = _early_stopping_patience // 10
_lr_scheduler: str = 'ReduceLROnPlateau'
_lr_scheduler_kwargs: Dict[str, Any] = {
    'factor': 0.2,
    'patience': _early_stopping_patience // 2,
    'threshold': 1e-3,
    'cooldown': 0,
}
# logging configurations
_num_trn_logs: int = 10

# plot and other documentation/visualization configurations
_emb_plot: bool = ((not _nni_search) and (not _dry_run))

# # adjust configurations if it's a dry run
# if _dry_run:
#     __dry_run_factor: int = 10
#     __warning_msg = \
#         f'Reducing the data length (sequences and paddings), ' \
#         f'the number of training and validation batches, and ' \
#         f'the number of training and early stopping patience epochs ' \
#         f'by a factor of {__dry_run_factor} for dry run ... '
#     __LOGGER.warning(__warning_msg)
#
#     _seq_len //= __dry_run_factor
#     _max_num_paddings = (_max_num_paddings // __dry_run_factor) \
#         if isinstance(_max_num_paddings, int) else _max_num_paddings
#     _max_num_trn_batches_per_epoch //= __dry_run_factor
#     _max_num_vld_batches_per_epoch //= __dry_run_factor
#     _max_num_tst_batches //= __dry_run_factor
#     _max_num_epochs //= __dry_run_factor
#     _early_stopping_patience //= __dry_run_factor


# dictionary that maps names of each configuration to their object
# e.g. 'experiment_name': _experiment_name, 'random_state': _random_state, etc.
# note that the local configuration variable names must start with '_', but
# the underscores are stripped away in the CONFIG dictionary
config: MappingProxyType = MappingProxyType({
    variable_name[1:]: variable
    for variable_name, variable in locals().items() if
    variable_name.startswith('_') and not variable_name.startswith('__')
})
