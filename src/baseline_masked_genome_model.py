"""
File Name:          baseline_masked_genome_model.py
Project:            bioseq-learning

File Description:

    Baseline model for masked genome prediction. References:
    https://github.com/pytorch/examples/blob/master/word_language_model/main.py
    https://github.com/bentrevett/pytorch-seq2seq/blob/master/6%20-%20Attention%20is%20All%20You%20Need.ipynb

"""
import copy
import time
import logging

import numpy as np
import torch
import torch.onnx
import torch.nn as nn
from torch.utils.data import DataLoader

from src import E_COLI_GENOME_PARENT_DIR_PATH
from src.datasets import split_genome_dir_paths, GenomeDataset
from src.datasets.genome_dataset import \
    PADDING_INDEX, NUCLEOTIDE_CHAR_INDEX_DICT, INDEX_NUCLEOTIDE_CHAR_DICT
from src.modules import TransformerEncoderModel
from src.optimization import get_torch_optimizer, get_torch_lr_scheduler
from src.utilities import set_random_seed, get_computation_devices
from src.configs.baseline_masked_genome_model_config import config


_LOGGER = logging.getLogger(__name__)

set_random_seed(
    random_seed=config['random_seed'],
    deterministic_cudnn_flag=config['deterministic_cudnn_flag'],
)
# TODO: multi-GPU model ...
device = get_computation_devices(
    preferred_gpu_list=config['preferred_gpu_list'],
    multi_gpu_flag=config['multi_gpu_flag'],
)[0]
# specify the cuda devices explicitly so that apex won't use cuda: 0
torch.cuda.set_device(device)

nvidia_amp_opt: bool = config['nvidia_amp_opt']
if config['nvidia_amp_opt']:
    try:
        from apex import amp
    except ImportError:
        _warning_msg = \
            f'Cannot import NVIDIA-apex. ' \
            f'Ignoring mixed-precision training/inference ... '
        _LOGGER.warning(_warning_msg)
        nvidia_amp_opt = False

# print out the configurations
print('=' * 80)
print('Configurations:')
print('-' * 80)
_config_key_max_len: int = max([len(_k) for _k in config.keys()])
for _key, _value in config.items():
    print(f'{_key:{_config_key_max_len + 4}s}: {_value}')
print('=' * 80)


###############################################################################
# Load data
###############################################################################

trn_genome_dir_paths, vld_genome_dir_paths, tst_genome_dir_paths = \
    split_genome_dir_paths(
        genome_parent_dir_path=E_COLI_GENOME_PARENT_DIR_PATH,
        vld_ratio=config['vld_ratio'],
        tst_ratio=config['tst_ratio']
    )
# if experimental (dry-run) indicator is set to True, then use only one
# genome, and setting validation/test sets to the same as the training one
trn_genome_dir_paths = \
    trn_genome_dir_paths[0:1] \
    if config['dry_run'] else trn_genome_dir_paths
masked_genome_dataset_kwargs = {
    'seq_len': config['seq_len'],
    # 'num_masks': config['num_masks'],
    'max_num_paddings': config['max_num_paddings'],
}
print(f'Generating training dataset from '
      f'{len(trn_genome_dir_paths)} genomes ...')
trn_dataset = GenomeDataset(
    trn_genome_dir_paths,
    **masked_genome_dataset_kwargs,
)
if config['dry_run']:
    vld_dataset = trn_dataset
    tst_dataset = trn_dataset
else:
    print(f'Generating validation dataset from '
          f'{len(vld_genome_dir_paths)} genomes ...')
    vld_dataset = GenomeDataset(
        vld_genome_dir_paths,
        **masked_genome_dataset_kwargs,
    )
    print(f'Generating testing dataset from '
          f'{len(tst_genome_dir_paths)} genomes ...')
    tst_dataset = GenomeDataset(
        tst_genome_dir_paths,
        **masked_genome_dataset_kwargs,
    )

dataloader_kwargs = {
    'batch_size': config['dataloader_batch_size'],
    'shuffle': True,
    'num_workers': config['dataloader_num_workers'],
    'pin_memory': (device.type == 'cuda'),
    # layer normalization requires that the input tensor has the same size,
    # and therefore requires dropping the last batch that might often be of
    # some different shape
    'drop_last': config['xfmr_enc_norm'],
}
trn_dataloader = DataLoader(trn_dataset, **dataloader_kwargs)
vld_dataloader = DataLoader(vld_dataset, **dataloader_kwargs)
tst_dataloader = DataLoader(tst_dataset, **dataloader_kwargs)

# the number of masks in each padded sequence
num_masks: int = int(np.round(config['num_masks'] * config['seq_len'])) \
    if isinstance(config['num_masks'], float) else config['num_masks']
# the number of tokens in the nucleotide char: index dict
num_tokens: int = len(NUCLEOTIDE_CHAR_INDEX_DICT)


###############################################################################
# Build the model
###############################################################################

model = TransformerEncoderModel(
    # number of tokens shall not include the padding token
    num_tokens=num_tokens,
    padding_index=PADDING_INDEX,
    seq_len=config['seq_len'],
    batch_size=config['dataloader_batch_size'],
    emb_dim=config['emb_dim'],
    pos_enc=config['pos_enc'],
    pos_enc_dropout=config['pos_enc_dropout'],
    pos_enc_emb_scale=config['pos_enc_emb_scale'],
    xfmr_enc_layer_num_attn_heads=config['xfmr_enc_layer_num_attn_heads'],
    xfmr_enc_layer_feedforward_dim=config['xfmr_enc_layer_feedforward_dim'],
    xfmr_enc_layer_activation=config['xfmr_enc_layer_activation'],
    xfmr_enc_layer_dropout=config['xfmr_enc_layer_dropout'],
    xfmr_enc_num_layers=config['xfmr_enc_num_layers'],
    xfmr_enc_norm=config['xfmr_enc_norm'],
).to(device)
criterion = nn.CrossEntropyLoss(ignore_index=PADDING_INDEX)
optimizer = get_torch_optimizer(
    optimizer=config['optimizer'],
    parameters=model.parameters(),
    optimizer_kwargs=config['optimizer_kwargs'],
)
if nvidia_amp_opt:
    print(f'Using Nvida-apex mixed-precision model. Configurations: ')
    model, optimizer = amp.initialize(
        model, optimizer, opt_level=config['nvidia_amp_opt_level'])

lr_scheduler = get_torch_lr_scheduler(
    lr_scheduler=config['lr_scheduler'],
    optimizer=optimizer,
    lr_scheduler_kwargs=config['lr_scheduler_kwargs'],
)


###############################################################################
# Training code
###############################################################################

# intervals between batches for training log
_num_trn_batches: int = \
    min(len(trn_dataloader), config['max_num_trn_batches_per_epoch'])
_trn_log_interval: int = int(_num_trn_batches / config['num_trn_logs'])


def train(cur_epoch: int):

    model.train()
    _total_loss = 0.
    _start_time = time.time()

    for _batch_index, _batch in enumerate(trn_dataloader):

        # stop this epoch if there has been too many batches already
        if _batch_index >= config['max_num_trn_batches_per_epoch']:
            break

        _indexed_seqs = _batch[0].transpose(0, 1).to(device)
        _key_padding_mask = _batch[1].to(device)

        # model.zero_grad()
        optimizer.zero_grad()

        _output, _ = \
            model(_indexed_seqs, num_masks, _key_padding_mask)
        _output = _output.view(-1, num_tokens)
        _target = _indexed_seqs.contiguous().view(-1)
        _loss = criterion(input=_output, target=_target)

        _loss.backward()
        _total_loss += _loss.item()
        optimizer.step()

        if (_batch_index + 1) % _trn_log_interval == 0:

            _avg_batch_loss: float = _total_loss / _trn_log_interval
            _avg_batch_time_in_ms: float = \
                (time.time() - _start_time) * 1000 / _trn_log_interval
            print(
                f'| epoch {cur_epoch:3d} '
                f'| {(_batch_index + 1):6d} / {_num_trn_batches:<d} batches '
                f'| learning rate {optimizer.param_groups[0]["lr"]:1.2E} '
                f'| loss {_avg_batch_loss:5.4f} '
                f'| {_avg_batch_time_in_ms:>4.0f} ms/batch |'
            )
            _total_loss = 0.
            _start_time = time.time()


def evaluate(_dataloader, test=False):

    model.eval()
    _total_loss = 0.
    _num_total_predictions = 0
    _num_correct_predictions = 0

    with torch.no_grad():
        for _batch_index, _batch in enumerate(trn_dataloader):

            # stop this epoch if there has been too many batches already
            # this is only applicable for validation
            if (not test) and \
                    _batch_index >= config['max_num_vld_batches_per_epoch']:
                break

            _indexed_seqs = _batch[0].transpose(0, 1).to(device)
            _key_padding_mask = _batch[1].to(device)

            _output, _mask = \
                model(_indexed_seqs, num_masks, _key_padding_mask)
            _output = _output.view(-1, num_tokens)
            _target = _indexed_seqs.contiguous().view(-1)
            _loss = criterion(input=_output, target=_target)

            # collect the metrics: loss and accuracy
            _total_loss += _loss.item()

            _mask = _mask.bool().repeat(config['dataloader_batch_size'])
            _, _prediction = torch.max(_output, 1)

            # TODO: turn this code segment into a printing function
            # __mask = _mask[:config['seq_len']].view(-1).tolist()
            # __target = _target[:config['seq_len']].view(-1).tolist()
            # __prediction = _prediction[:config['seq_len']].view(-1).tolist()
            #
            # # print(__mask)
            #
            # for _i in range(len(__target)):
            #     print(INDEX_NUCLEOTIDE_CHAR_DICT[__target[_i]], end='')
            # print('')
            # for _i in range(len(__mask)):
            #     if __mask[_i]:
            #         print('|', end='')
            #     else:
            #         print(' ', end='')
            # print('')
            # for _i in range(len(__prediction)):
            #     print(INDEX_NUCLEOTIDE_CHAR_DICT[__prediction[_i]], end='')
            # print('')

            _masked_target = _target[_mask]
            _masked_prediction = _prediction[_mask]
            _no_padding_masks = _masked_target != PADDING_INDEX

            _num_total_predictions += _no_padding_masks.sum().item()
            _num_correct_predictions += (
                    _no_padding_masks &
                    (_masked_prediction == _masked_target)
            ).sum().item()

    _num_batches: int = min(
        len(_dataloader),
        config['max_num_vld_batches_per_epoch'],
    )
    _loss = _total_loss / _num_batches
    _acc = _num_correct_predictions / _num_total_predictions
    return _loss, _acc


# train the model over the epochs and evaluate on the validation set
best_vld_loss, best_vld_acc, best_epoch, best_model = \
    float('inf'), 0., 0, None
print('=' * 80)
try:
    for epoch in range(1, config['max_num_epochs']):

        epoch_start_time = time.time()
        train(epoch)
        epoch_vld_loss, epoch_vld_acc = evaluate(vld_dataloader)
        lr_scheduler.step()
        epoch_time_in_sec = time.time() - epoch_start_time

        print('-' * 80)
        print(
            f'| end of epoch {epoch:3d} '
            f'| time {epoch_time_in_sec:>5.0f} s '
            f'| validation loss {epoch_vld_loss:5.4f} '
            f'| validation accuracy {(epoch_vld_acc * 100):3.3f}% '
        )
        print('-' * 80)

        # if epoch_vld_loss < best_vld_loss:
        if epoch_vld_acc > best_vld_acc:
            best_vld_loss = epoch_vld_loss
            best_vld_acc = epoch_vld_acc
            best_epoch = epoch
            best_model = copy.deepcopy(model)
        elif epoch - best_epoch >= config['early_stopping_patience']:
            print('exiting from training early for early stopping ... ')

except KeyboardInterrupt:
    print('exiting from training early for KeyboardInterrupt ... ')


# evaluate the model on the test set
if best_model:
    model = best_model
    tst_loss, tst_acc = evaluate(tst_dataloader, test=True)
    print('=' * 80)
    print(
        f'| End of training '
        f'| test loss {tst_loss:5.4f} '
        f'| test accuracy {(tst_acc * 100):3.2f}% '
    )
    print('=' * 80)
