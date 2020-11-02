"""
File Name:          baseline_masked_genome_model.py
Project:            bioseq-learning

File Description:

    Baseline model for masked genome prediction. References:
    https://github.com/pytorch/examples/blob/master/word_language_model/main.py
    https://github.com/bentrevett/pytorch-seq2seq/blob/master/6%20-%20Attention%20is%20All%20You%20Need.ipynb

"""
import os
import sys
import math
import time
import logging
import traceback
from typing import Any, Dict
from types import MappingProxyType

import numpy as np
import torch
import torch.onnx
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.optim.lr_scheduler import ReduceLROnPlateau
from pytorch_model_summary import summary

from src import E_COLI_GENOME_PARENT_DIR_PATH
from src.datasets.split_genome_dir_paths import split_genome_dir_paths
from src.datasets.genome_dataset import \
    PADDING_INDEX, NUCLEOTIDE_CHAR_INDEX_DICT, \
    GenomeDataset, GenomeIterDataset
from src.modules.transformer_encoder_model import \
    Masker, Seq2Kmer, Kmer2Seq, get_transformer_encoder_model
from src.optimization import get_torch_optimizer, get_torch_lr_scheduler
from src.utilities import set_random_seed, get_computation_devices
from src.utilities.merge_nni_config import merge_nni_config
from src.configs.baseline_masked_genome_model_config import config as \
    default_config


_LOGGER = logging.getLogger(__name__)


# configure the trial parameters with given configuration and NNI settings
if default_config['nni_search']:
    import nni

    # merge the configured hyper-parameters and the next tuning parameters
    # notes on the NNI config v1.8.0:
    # - all the values are numeric or strings (of choices)
    # - nested configurations are stored in list instead of dict
    _nni_config: Dict[str, Any] = nni.get_next_parameter()
    assert len(_nni_config) > 0

    config: MappingProxyType = \
        merge_nni_config(default_config, _nni_config)

    # sanity check the configurations to avoid failed runs after building
    # the datasets, models, and the optimizer etc.
    assert config['emb_dim'] % 2 == 0
    assert config['emb_dim'] % config['xfmr_enc_layer_num_attn_heads'] == 0

    # NNI will automatically configure the GPU(s) assigned to this trial
    # note that torch.cuda.current_device() will always give you 0
    devices = get_computation_devices(
        preferred_gpu_list='all',
        multi_gpu_flag=config['multi_gpu_flag'],
    )
    nni_search: bool = True
    checkpoint_dir_path: str = os.path.join(
        config['model_directory'],
        f'nni/{nni.get_experiment_id()}'
    )
    os.makedirs(checkpoint_dir_path, exist_ok=True)
    checkpoint_path: str = os.path.join(
        checkpoint_dir_path,
        f'{nni.get_trial_id()}.pt',
    )
else:
    config: MappingProxyType = default_config
    devices = get_computation_devices(
        preferred_gpu_list=config['preferred_gpu_list'],
        multi_gpu_flag=config['multi_gpu_flag'],
    )
    nni_search: bool = False
    checkpoint_path: str = os.path.join(
        config['model_directory'], f'temporary.pt')

set_random_seed(
    random_seed=config['random_seed'],
    deterministic_cudnn_flag=config['deterministic_cudnn_flag'],
)
if devices[0].type == 'cuda':
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
        # note that if device number is not explicitly set, then the following
        # line won't work ... which essentially means tht apex and nni cannot
        # work well together, in the way that apex might always use a small
        # chunk of memory on cuda:0 even if the device is not cuda:0
        except ValueError:
            _warning_msg = \
                f'Failed tp specify cuda device; and Apex might take up ' \
                f'a small chunk of memory in cuda:0 during the run.'
else:
    nvidia_amp_opt: bool = False


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

# the number of masks in each padded sequence
num_masks: int = int(np.round(config['num_masks'] * config['seq_len'])) \
    if isinstance(config['num_masks'], float) else config['num_masks']
# the number of tokens in the nucleotide char: index dict
num_tokens: int = len(NUCLEOTIDE_CHAR_INDEX_DICT)
num_kmer_tokens: int = (num_tokens ** config['kmer_len'])

trn_genome_dir_paths, vld_genome_dir_paths, tst_genome_dir_paths = \
    split_genome_dir_paths(
        genome_parent_dir_path=E_COLI_GENOME_PARENT_DIR_PATH,
        vld_ratio=config['vld_ratio'],
        tst_ratio=config['tst_ratio']
    )
# if experimental (dry-run) indicator is set to True, then use only one
# genome for training, validation and testing separately
if config['dry_run']:
    # could pick genome 1033813.3 for one-contig genome
    # use the same training and validation set to see if the model works
    trn_genome_dir_paths = trn_genome_dir_paths[:10]
    vld_genome_dir_paths = trn_genome_dir_paths
    tst_genome_dir_paths = trn_genome_dir_paths

_max_num_paddings: int = config['max_num_paddings'] \
    if isinstance(config['max_num_paddings'], int) else \
    int(config['max_num_paddings'] * config['seq_len'])
masked_genome_dataset_kwargs = {
    'seq_len': config['seq_len'],
    'max_num_paddings': _max_num_paddings,
    'padding_mask': config['xfmr_padding_mask'],
    'dry_run': config['dry_run'],
    'dry_run_contig_len': config['dry_run_contig_len'],
    'dry_run_num_contigs': config['dry_run_num_contigs'],
    'dry_run_num_ins_per_contig': config['dry_run_num_ins_per_contig'],
}
print(f'Generating training dataset from '
      f'{len(trn_genome_dir_paths)} genomes ...')
trn_iter_dataset = GenomeIterDataset(
    trn_genome_dir_paths,
    strict_iteration=False,
    **masked_genome_dataset_kwargs,
)
print(f'Generating validation dataset from '
      f'{len(vld_genome_dir_paths)} genomes ...')
vld_iter_dataset = GenomeIterDataset(
    vld_genome_dir_paths,
    strict_iteration=False,
    **masked_genome_dataset_kwargs,
) if (vld_genome_dir_paths != trn_genome_dir_paths) else trn_iter_dataset
print(f'Generating testing dataset from '
      f'{len(tst_genome_dir_paths)} genomes ...')
tst_dataset = GenomeDataset(
    tst_genome_dir_paths,
    **masked_genome_dataset_kwargs,
)


# modified collection function so that sequence shape is compatible with
# transformer (sequence_length, batch_size)
# def __collate_fn(__batch):
#     __batch = default_collate(__batch)
#     __indexed_seqs = __batch[0].transpose(0, 1).contiguous()
#     __padding_mask = __batch[1]
#     return __indexed_seqs, __padding_mask


dataloader_kwargs = {
    'batch_size': config['dataloader_batch_size'],
    'num_workers': config['dataloader_num_workers'],
    # 'pin_memory': (devices[0].type == 'cuda'),
    # shuffle should be set to False for training adn validation:
    # (0) shuffling all the genome sequences takes extremely long and a good
    #     chunk of memory (> 150 G)
    # (1) training and validation sets are both iterable datasets,
    #     which requires no shuffling in the dataloader level
    # 'shuffle': False,
    #
    # set 'drop_last' to True all the time to avoid trouble caused by
    # irregular batch size for normalization and GPipe chunking
    'drop_last': True,
    # 'collate_fn': __collate_fn,
}

trn_dataloader = DataLoader(trn_iter_dataset, **dataloader_kwargs)
vld_dataloader = DataLoader(vld_iter_dataset, **dataloader_kwargs)
tst_dataloader = DataLoader(tst_dataset, shuffle=True, **dataloader_kwargs)


###############################################################################
# Build the model
###############################################################################

masker = Masker(
    seq_len=config['seq_len'],
    num_masks=num_masks,
    mask_index=PADDING_INDEX,
    attn_mask=config['xfmr_attn_mask'],
)
seq2kmer = Seq2Kmer(
    k=config['kmer_len'],
    seq_len=config['seq_len'],
    num_tokens=num_tokens,
    attn_mask=config['xfmr_attn_mask'],
    padding_mask=config['xfmr_padding_mask'],
)
model = get_transformer_encoder_model(
    # number of tokens shall not include the padding token
    num_tokens=num_kmer_tokens,
    padding_index=PADDING_INDEX,
    seq_len=(config['seq_len'] - config['kmer_len'] + 1),
    emb_dim=config['emb_dim'],
    pos_enc=config['pos_enc'],
    pos_enc_dropout=config['pos_enc_dropout'],
    pos_enc_emb_scale=math.sqrt(config['emb_dim']),
    xfmr_enc_layer_num_attn_heads=config['xfmr_enc_layer_num_attn_heads'],
    xfmr_enc_layer_feedforward_dim=config['xfmr_enc_layer_feedforward_dim'],
    xfmr_enc_layer_activation=config['xfmr_enc_layer_activation'],
    xfmr_enc_layer_dropout=config['xfmr_enc_layer_dropout'],
    xfmr_enc_layer_norm=config['xfmr_enc_layer_norm'],
    xfmr_enc_num_layers=config['xfmr_enc_num_layers'],
    xfmr_attn_mask=config['xfmr_attn_mask'],
    xfmr_padding_mask=config['xfmr_padding_mask'],
)
kmer2seq = Kmer2Seq(
    k=config['kmer_len'],
    seq_len=config['seq_len'],
    num_tokens=num_tokens,
)

# print model with sampled input
_batch = next(iter(vld_dataloader))
_masked_indexed_seqs, _attn_mask = masker(_batch[0])
_input_sample = seq2kmer((_masked_indexed_seqs, _attn_mask, _batch[1]))
print(summary(model, _input_sample, show_input=True))

if config['multi_gpu_flag']:
    from src.utilities.shard_module import shard_module
    # _batch = next(iter(vld_dataloader))
    # _masked_indexed_seqs, _attn_mask = masker(_batch[0])
    # _input_sample = seq2kmer((_masked_indexed_seqs, _attn_mask, _batch[1]))
    model, devices = shard_module(
        module=model,
        input_sample=_input_sample,
        devices=devices,
        # if attention masks are passing around between GPUs, it should not
        # be chunked in by the pipeline as it's not correlated to the batch
        # size. Should i simply use broadcast instead?
        num_chunks=1 if config['xfmr_attn_mask']
        else config['num_gpipe_chunks'],
    )
    masker = masker.to(devices[0])
    seq2kmer = seq2kmer.to(devices[0])
    kmer2seq = kmer2seq.to(devices[-1])
else:
    # if the multi-gpu flag is False/disabled
    masker = masker.to(devices[0])
    seq2kmer = seq2kmer.to(devices[0])
    model = model.to(devices[0])
    kmer2seq = kmer2seq.to(devices[0])

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
    config['max_num_trn_batches_per_epoch'] if config['dry_run'] else \
    min(len(trn_dataloader), config['max_num_trn_batches_per_epoch'])
_trn_log_interval: int = int(_num_trn_batches / config['num_trn_logs'])


def train(cur_epoch: int):

    model.train()
    _total_trn_loss = 0.
    _start_time = time.time()

    for _batch_index, _batch in enumerate(trn_dataloader):

        # stop this epoch if there has been too many batches already
        if _batch_index >= _num_trn_batches:
            break

        _indexed_seqs = _batch[0].to(devices[0], non_blocking=True)
        # using an empty tensor for padding mask if the padding mask is 
        # disabled, to save the memory and IO bandwidth
        _padding_mask = \
            _batch[1].to(devices[0], non_blocking=True) \
            if config['xfmr_padding_mask'] else torch.zeros(
                size=(config['dataloader_batch_size'], 0),
                dtype=torch.bool).to(devices[0], non_blocking=True)

        # model.zero_grad()
        optimizer.zero_grad()

        # update the mask for every batch
        # if not config['dry_run']:
        #     masker.update()
        masker.update()

        # TODO: switch some data processing from module to dataset
        # the model output has the shape of (batch_size, seq_len, num_tokens)
        _masked_indexed_seqs, _attn_mask = masker(_indexed_seqs)
        _indexed_kmer_seqs, _, _ = seq2kmer((_indexed_seqs, None, None))
        _masked_indexed_kmer_seqs, _attn_kmer_mask, _padding_kmer_mask = \
            seq2kmer((_masked_indexed_seqs, _attn_mask, _padding_mask))
        _output = model((
            _masked_indexed_kmer_seqs,
            _attn_kmer_mask,
            _padding_kmer_mask,
        ))

        _trn_loss = criterion(
            input=_output.view(-1, num_kmer_tokens),
            target=_indexed_kmer_seqs.
                to(devices[-1], non_blocking=True).view(-1),
        )
        _trn_loss.backward()
        nn.utils.clip_grad_norm_(
            parameters=model.parameters(),
            max_norm=config['max_grad_norm'],
        )
        optimizer.step()

        _trn_loss: float = _trn_loss.item()
        _total_trn_loss += _trn_loss

        if (_batch_index + 1) % _trn_log_interval == 0:

            _trn_avg_loss: float = _total_trn_loss / _trn_log_interval
            _trn_avg_time_in_ms: float = \
                (time.time() - _start_time) * 1000 / _trn_log_interval
            print(
                f'| epoch {cur_epoch:3d} '
                f'| {(_batch_index + 1):6d} / {_num_trn_batches:<d} batches '
                f'| learning rate {optimizer.param_groups[0]["lr"]:1.2E} '
                f'| loss {_trn_avg_loss:7.4f} '
                f'| {_trn_avg_time_in_ms:>4.0f} ms/batch |'
            )

            _total_trn_loss = 0.
            _start_time = time.time()


def evaluate(_dataloader, test=False):

    model.eval()
    _total_loss = 0.
    _num_total_predictions = 0
    _num_correct_predictions = 0

    max_num_batches: int = \
        min(len(tst_dataloader), config['max_num_tst_batches']) if test \
        else min(len(vld_dataloader), config['max_num_vld_batches_per_epoch'])

    with torch.no_grad():
        for _batch_index, _batch in enumerate(_dataloader):

            # stop this epoch if there has been too many batches already
            if _batch_index >= max_num_batches:
                break

            _indexed_seqs = _batch[0].to(devices[0], non_blocking=True)
            _padding_mask = \
                _batch[1].to(devices[0], non_blocking=True) \
                if config['xfmr_padding_mask'] else torch.zeros(
                    size=(config['dataloader_batch_size'], 0),
                    dtype=torch.bool).to(devices[0], non_blocking=True)

            # if not config['dry_run']:
            #     masker.update()
            masker.update()

            _masked_indexed_seqs, _attn_mask = masker(_indexed_seqs)
            _indexed_kmer_seqs, _, _ = seq2kmer((_indexed_seqs, None, None))
            _masked_indexed_kmer_seqs, _attn_kmer_mask, _padding_kmer_mask = \
                seq2kmer((_masked_indexed_seqs, _attn_mask, _padding_mask))
            _output = model((
                _masked_indexed_kmer_seqs,
                _attn_kmer_mask,
                _padding_kmer_mask,
            ))

            _input = _masked_indexed_kmer_seqs.view(-1).to(
                devices[-1], non_blocking=True)
            _target = _indexed_kmer_seqs.view(-1).to(
                devices[-1], non_blocking=True)
            _output = _output.view(-1, num_kmer_tokens)

            # collect the loss
            _loss = criterion(input=_output, target=_target)
            _total_loss += _loss.item()

            # collect the accuracy
            _, _prediction = torch.max(_output, 1)

            # prediction masks over the whole batch is defined by where the
            # input and target sequences differ; in this case, the masked
            # padding is still a mask, therefore will not be counted
            # note that here is how to get the prediction masks over the
            # whole batch including the paddings:
            # _batch_mask = _src_mask.repeat(config['dataloader_batch_size'])
            # _batch_mask = _batch_mask.view(-1, config['seq_len'])
            # _batch_mask = _batch_mask.transpose(1, 0).contiguous()
            # _batch_mask = _batch_mask.view(-1)
            _batch_mask: torch.BoolTensor = (_input != _target)

            _masked_target = _target[_batch_mask]
            _masked_prediction = _prediction[_batch_mask]
            _num_total_predictions += _batch_mask.sum().item()
            _num_correct_predictions += \
                (_masked_prediction == _masked_target).sum().item()

            # # print the predictions and the targets
            # print_masked_genome_predictions(
            #     config['seq_len'],
            #     config['dataloader_batch_size'],
            #     _input, _target, _prediction,
            # )

    _num_batches: int = \
        min(len(_dataloader), config['max_num_tst_batches']) if test else \
        min(len(_dataloader), config['max_num_vld_batches_per_epoch'])
    _loss = _total_loss / _num_batches
    _acc = _num_correct_predictions / _num_total_predictions
    return _loss, _acc


if __name__ == '__main__':

    # train the model over the epochs and evaluate on the validation set
    best_vld_loss, best_vld_acc, best_epoch = float('inf'), 0., 0
    num_epochs: int = 0
    print('=' * 80)
    while True:
        try:
            for epoch in range(1, config['max_num_epochs']):
                num_epochs += 1
                epoch_start_time = time.time()
                train(epoch)
                epoch_vld_loss, epoch_vld_acc = evaluate(vld_dataloader)
                if nni_search:
                    nni.report_intermediate_result({
                        'default': epoch_vld_acc,
                        # 'vld_acc': epoch_vld_acc,
                        'vld_avg_loss': epoch_vld_loss,
                    })

                if isinstance(lr_scheduler, ReduceLROnPlateau):
                    lr_scheduler.step(metrics=epoch_vld_loss)
                else:
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
                # checkpoint model if it has the best performance so far
                if epoch_vld_acc > best_vld_acc:
                    best_vld_loss = epoch_vld_loss
                    best_vld_acc = epoch_vld_acc
                    best_epoch = epoch
                    # best_model = copy.deepcopy(model)
                    torch.save(
                        {
                            'epoch': epoch,
                            'model_state_dict': model.state_dict(),
                            'optimizer_state_dict': optimizer.state_dict(),
                        },
                        checkpoint_path,
                    )
                elif epoch - best_epoch >= config['early_stopping_patience']:
                    print('exiting from training early '
                          'for early stopping ... ')
                    break

                if config['dry_run'] and (epoch_vld_acc > 0.8):
                    print('exiting from dry run training early '
                          'for > 80% accuracy ... ')
                    break

                if math.isnan(epoch_vld_loss):
                    print('validation loss gets to NaN; exiting current '
                          'invalid trial ...')
                    # if nni_search:
                    #     nni.report_final_result({'default': 0})
                    sys.exit(0)
            break

        except RuntimeError as e:
            # CUDA memory or other errors from the training/evaluation process
            traceback.print_exc()
            # if nni_search:
            #     nni.report_final_result({'default': 0})
            sys.exit(1)

        except KeyboardInterrupt:
            print('exiting from training early for KeyboardInterrupt ... ')
            break

    # evaluate the model on the test set
    model.load_state_dict(torch.load(checkpoint_path)['model_state_dict'])
    tst_start_time = time.time()
    tst_loss, tst_acc = evaluate(tst_dataloader, test=True)
    tst_time_in_sec = time.time() - tst_start_time

    if nni_search:
        nni.report_final_result({
            'default': tst_acc,
            # 'tst_acc': tst_acc,
            'num_epochs': num_epochs,
            'tst_loss': tst_loss,
        })

    print('=' * 80)
    print(
        f'| end of training '
        f'| test time {tst_time_in_sec:>5.0f} s '
        f'| test loss {tst_loss:5.4f} '
        f'| test accuracy {(tst_acc * 100):3.2f}% '
    )
    print('=' * 80)
