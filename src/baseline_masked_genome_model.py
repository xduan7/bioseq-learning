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

from src import E_COLI_GENOME_PARENT_DIR_PATH
from src.datasets.split_genome_dir_paths import split_genome_dir_paths
from src.datasets.genome_dataset import \
    PADDING_INDEX, NUCLEOTIDE_CHAR_INDEX_DICT, \
    GenomeDataset, GenomeIterDataset
from src.modules.transformer_encoder_model import get_transformer_encoder_model
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

    # NNI will automatically configure the GPU(s) assigned to this trial
    # note that torch.cuda.current_device() will always give you 0
    devices = get_computation_devices(
        preferred_gpu_list='all',
        multi_gpu_flag=config['multi_gpu_flag'],
    )
    nni_search: bool = True
    model_checkpoint_path: str = os.path.join(
        config['model_directory'],
        f'nni-{nni.get_trial_id()}.pt',
    )
else:
    config: MappingProxyType = default_config
    devices = get_computation_devices(
        preferred_gpu_list=config['preferred_gpu_list'],
        multi_gpu_flag=config['multi_gpu_flag'],
    )
    nni_search: bool = False
    model_checkpoint_path: str = os.path.join(
        config['model_directory'], f'temporary.pt')

# configure computation devices here ...
# TODO: multi-GPU model ...


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
    trn_genome_dir_paths = trn_genome_dir_paths[0:10]
    vld_genome_dir_paths = vld_genome_dir_paths[0:10]
    tst_genome_dir_paths = tst_genome_dir_paths[0:1]

masked_genome_dataset_kwargs = {
    'seq_len': config['seq_len'],
    'max_num_paddings': config['max_num_paddings'],
}
print(f'Generating training dataset from '
      f'{len(trn_genome_dir_paths)} genomes ...')
trn_iter_dataset = GenomeIterDataset(
    trn_genome_dir_paths,
    **masked_genome_dataset_kwargs,
)
print(f'Generating validation dataset from '
      f'{len(vld_genome_dir_paths)} genomes ...')
vld_iter_dataset = GenomeIterDataset(
    vld_genome_dir_paths,
    **masked_genome_dataset_kwargs,
)
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
    'pin_memory': (devices[0].type == 'cuda'),
    # shuffle should be set to False for training adn validation:
    # (0) shuffling all the genome sequences takes extremely long and a good
    #     chunk of memory (> 150 G)
    # (1) training and validation sets are both iterable datasets,
    #     which requires no shuffling in the dataloader level
    # 'shuffle': False,
    #
    # layer normalization requires that the input tensor has the same size
    # therefore requires dropping the last batch of data, which is often of
    # different shapes compared to previous batches
    'drop_last': config['xfmr_enc_layer_norm'],
    # 'collate_fn': __collate_fn,
}

trn_dataloader = DataLoader(trn_iter_dataset, **dataloader_kwargs)
vld_dataloader = DataLoader(vld_iter_dataset, **dataloader_kwargs)
tst_dataloader = DataLoader(tst_dataset, shuffle=True, **dataloader_kwargs)



###############################################################################
# Build the model
###############################################################################

model = get_transformer_encoder_model(
    # number of tokens shall not include the padding token
    num_tokens=num_tokens,
    num_masks=num_masks,
    padding_index=PADDING_INDEX,
    mask_index=PADDING_INDEX,
    seq_len=config['seq_len'],
    emb_dim=config['emb_dim'],
    pos_enc=config['pos_enc'],
    pos_enc_dropout=config['pos_enc_dropout'],
    pos_enc_emb_scale=config['pos_enc_emb_scale'],
    xfmr_enc_layer_num_attn_heads=config['xfmr_enc_layer_num_attn_heads'],
    xfmr_enc_layer_feedforward_dim=config['xfmr_enc_layer_feedforward_dim'],
    xfmr_enc_layer_activation=config['xfmr_enc_layer_activation'],
    xfmr_enc_layer_dropout=config['xfmr_enc_layer_dropout'],
    xfmr_enc_layer_norm=config['xfmr_enc_layer_norm'],
    xfmr_enc_num_layers=config['xfmr_enc_num_layers'],
)
if config['multi_gpu_flag']:

    from torchgpipe import GPipe
    from torchgpipe.balance.blockpartition import solve
    from src.utilities.get_module_summary import get_module_summary

    # summarize the model with a forward pass
    # validation dataloader is faster compared to training (bigger) or test
    # (strictly mapping) datasets
    _batch = next(iter(vld_dataloader))
    _model_summary_ordered_dict, _model_summary_str = \
        get_module_summary(model, tuple(_batch))

    # get a dict that maps sequential layer name -> number of parameters
    _model_summary_list = list(_model_summary_ordered_dict.items())
    _summary_index = 0
    _layer_num_params_dict: Dict[str, int] = {}
    for _layer in model:

        _layer_num_params = 0
        _layer_class = str(_layer.__class__).split('.')[-1].split('\'')[0]

        while not _model_summary_list[_summary_index][0]\
                .startswith(_layer_class):
            _layer_num_params += \
                _model_summary_list[_summary_index][1]['num_params']
            _summary_index += 1

        _summary_index += 1

        __i: int = 0
        for __k in _layer_num_params_dict.keys():
            if __k.startswith(_layer_class):
                __i += 1
        _layer_name = _layer_class + f'-{__i}'
        _layer_num_params_dict[_layer_name] = _layer_num_params

    # TODO: check how many devices are sufficient enough
    # reduce the 'devices' list to the minimal number of devices
    # torch.cuda.get_device_properties(0).total_memory

    # solve for balanced distribution over all the devices
    _balance = [len(__b) for __b in solve(
        list(_layer_num_params_dict.values()),
        partitions=len(devices),
    )]

    #
    model = GPipe(
        module=model,
        balance=_balance,
        devices=devices,
        chunks=1,
        checkpoint='never',
    )

else:
    model = model.to(devices[0])

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
    _total_trn_loss = 0.
    _start_time = time.time()

    for _batch_index, _batch in enumerate(trn_dataloader):

        # stop this epoch if there has been too many batches already
        if _batch_index >= config['max_num_trn_batches_per_epoch']:
            break

        _indexed_seqs = _batch[0].to(devices[0], non_blocking=True)
        _padding_mask = _batch[1].to(devices[0], non_blocking=True)

        # model.zero_grad()
        optimizer.zero_grad()

        # update the mask for every batch
        model[0].update()

        # the model output has the shape of (batch_size, seq_len, num_tokens)
        _output = model((_indexed_seqs, _padding_mask))

        _trn_loss = criterion(
            input=_output.view(-1, num_tokens),
            target=_indexed_seqs.to(devices[-1], non_blocking=True).view(-1),
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
                f'| loss {_trn_avg_loss:5.4f} '
                f'| {_trn_avg_time_in_ms:>4.0f} ms/batch |'
            )

            _total_trn_loss = 0.
            _start_time = time.time()


def evaluate(_dataloader, test=False):

    model.eval()
    _total_loss = 0.
    _num_total_predictions = 0
    _num_correct_predictions = 0

    max_num_batches: int = config['max_num_tst_batches'] if test \
        else config['max_num_vld_batches_per_epoch']

    with torch.no_grad():
        for _batch_index, _batch in enumerate(_dataloader):

            # stop this epoch if there has been too many batches already
            if _batch_index >= max_num_batches:
                break

            _indexed_seqs = _batch[0].to(devices[0], non_blocking=True)
            _padding_mask = _batch[1].to(devices[0], non_blocking=True)

            model[0].update()
            _output = model((_indexed_seqs, _padding_mask))

            _input = model[0].curr_masked_indexed_seqs.view(-1).to(
                devices[-1], non_blocking=True)
            _target = _indexed_seqs.view(-1).to(
                devices[-1], non_blocking=True)
            _output = _output.view(-1, num_tokens)

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
    print('=' * 80)
    while True:
        try:
            for epoch in range(1, config['max_num_epochs']):

                epoch_start_time = time.time()
                train(epoch)
                epoch_vld_loss, epoch_vld_acc = evaluate(vld_dataloader)
                if nni_search:
                    nni.report_intermediate_result({
                        'default': epoch_vld_acc,
                        # 'vld_acc': epoch_vld_acc,
                        'vld_avg_loss': epoch_vld_loss,
                    })

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
                    # best_model = copy.deepcopy(model)
                    torch.save(model.state_dict(), model_checkpoint_path)
                elif epoch - best_epoch >= config['early_stopping_patience']:
                    print('exiting from training early for early stopping ... ')
                    break

                if math.isnan(epoch_vld_loss):
                    print('validation loss gets to NaN; exiting current '
                          'invalid trial ...')
                    # if nni_search:
                    #     nni.report_final_result({'default': 0})
                    sys.exit()
            break

        except RuntimeError as e:
            # CUDA memory or other errors from the training/evaluation process
            traceback.print_exc()
            # if nni_search:
            #     nni.report_final_result({'default': 0})
            sys.exit()

        except KeyboardInterrupt:
            print('exiting from training early for KeyboardInterrupt ... ')
            break

    # evaluate the model on the test set
    model.load_state_dict(torch.load(model_checkpoint_path))
    tst_start_time = time.time()
    tst_loss, tst_acc = evaluate(tst_dataloader, test=True)
    tst_time_in_sec = time.time() - tst_start_time

    if nni_search:
        nni.report_final_result({
            'default': tst_acc,
            # 'tst_acc': tst_acc,
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
