"""
File Name:          baseline_masked_genome_model.py
Project:            bioseq-learning

File Description:

    Baseline model for masked genome prediction. References:
    https://github.com/pytorch/examples/blob/master/word_language_model/main.py
    https://github.com/bentrevett/pytorch-seq2seq/blob/master/6%20-%20Attention%20is%20All%20You%20Need.ipynb

"""
import sys
import copy
import math
import time
import logging
import traceback

import numpy as np
import torch
import torch.onnx
import torch.nn as nn
from torch.utils.data import DataLoader
from torch.utils.data.dataloader import default_collate

from src import E_COLI_GENOME_PARENT_DIR_PATH
from src.datasets import \
    split_genome_dir_paths, print_masked_genome_predictions, \
    SequenceMask, GenomeDataset, GenomeIterDataset
from src.datasets.genome_dataset import \
    PADDING_INDEX, NUCLEOTIDE_CHAR_INDEX_DICT
from src.modules import TransformerEncoderModel
from src.optimization import get_torch_optimizer, get_torch_lr_scheduler
from src.utilities import set_random_seed, get_computation_devices
from src.configs.baseline_masked_genome_model_config import config as \
    default_config


_LOGGER = logging.getLogger(__name__)


# configure the trial parameters with given configuration and NNI settings
if default_config['nni_search']:
    import nni
    # merge the configured hyper-parameters and the next tuning parameters
    # TODO: better merge function
    #  (1) matches dtype
    #  (2) clarify which one is base
    config = {**dict(default_config), **nni.get_next_parameter()}

    # needs to do since NNI config cannot accept anything else other than
    # numbers or strings, and therefore if there is hyper-param of any other
    # type, one must manually convert the type
    config['xfmr_enc_norm']: bool = config['xfmr_enc_norm'] if \
        isinstance(config['xfmr_enc_norm'], bool) else \
        (config['xfmr_enc_norm'].lower() == 'true')

    # NNI will automatically configure the GPU(s) assigned to this trial
    # device = get_computation_devices(
    #     preferred_gpu_list=config['preferred_gpu_list'],
    #     multi_gpu_flag=config['multi_gpu_flag'],
    # )[0]
    device = torch.device(torch.cuda.current_device())
    nni_search: bool = True
else:
    config = default_config
    # TODO: multi-GPU model ...
    device = get_computation_devices(
        preferred_gpu_list=config['preferred_gpu_list'],
        multi_gpu_flag=config['multi_gpu_flag'],
    )[0]
    nni_search: bool = False

set_random_seed(
    random_seed=config['random_seed'],
    deterministic_cudnn_flag=config['deterministic_cudnn_flag'],
)
if device.type == 'cuda':
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
def __collate_fn(__batch):
    __batch = default_collate(__batch)
    __indexed_seqs = __batch[0].transpose(0, 1).contiguous()
    __padding_mask = __batch[1]
    return __indexed_seqs, __padding_mask


dataloader_kwargs = {
    'batch_size': config['dataloader_batch_size'],
    'num_workers': config['dataloader_num_workers'],
    'pin_memory': (device.type == 'cuda'),
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
    'drop_last': config['xfmr_enc_norm'],
    'collate_fn': __collate_fn,
}

trn_dataloader = DataLoader(trn_iter_dataset, **dataloader_kwargs)
vld_dataloader = DataLoader(vld_iter_dataset, **dataloader_kwargs)
tst_dataloader = DataLoader(tst_dataset, shuffle=True, **dataloader_kwargs)

# create a SequenceMask which generates the sequence and attention masks
seq_mask = SequenceMask(seq_len=config['seq_len'], num_masks=num_masks)


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
    _total_trn_loss = 0.
    _start_time = time.time()

    for _batch_index, _batch in enumerate(trn_dataloader):

        # stop this epoch if there has been too many batches already
        if _batch_index >= config['max_num_trn_batches_per_epoch']:
            break

        _indexed_seqs = _batch[0].to(device)
        _padding_mask = _batch[1].to(device)

        _attn_mask, _src_mask = seq_mask.update()
        _attn_mask = _attn_mask.to(device)
        _src_mask = _src_mask.to(device)

        _masked_indexed_seqs = _indexed_seqs.clone()
        _masked_indexed_seqs[_src_mask, :] = 0

        # model.zero_grad()
        optimizer.zero_grad()

        # the model output has the shape of (seq_len, batch_size, num_tokens)
        _output = model(
            src=_masked_indexed_seqs,
            attn_mask=_attn_mask,
            padding_mask=_padding_mask,
        )
        _trn_loss = criterion(
            input=_output.view(-1, num_tokens),
            target=_indexed_seqs.view(-1),
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

            _indexed_seqs = _batch[0].to(device)
            _padding_mask = _batch[1].to(device)

            _attn_mask, _src_mask = seq_mask.update()
            _attn_mask = _attn_mask.to(device)
            _src_mask = _src_mask.to(device)

            _masked_indexed_seqs = _indexed_seqs.clone()
            _masked_indexed_seqs[_src_mask, :] = 0

            _output = model(
                src=_masked_indexed_seqs,
                attn_mask=_attn_mask,
                padding_mask=_padding_mask,
            )
            _input = _masked_indexed_seqs.view(-1)
            _output = _output.view(-1, num_tokens)
            _target = _indexed_seqs.view(-1)

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

    _num_batches: int = len(_dataloader) if test else \
        min(len(_dataloader), config['max_num_vld_batches_per_epoch'])
    _loss = _total_loss / _num_batches
    _acc = _num_correct_predictions / _num_total_predictions
    return _loss, _acc


# train the model over the epochs and evaluate on the validation set
best_vld_loss, best_vld_acc, best_epoch, best_model = \
    float('inf'), 0., 0, None
print('=' * 80)
while True:
    try:
        for epoch in range(1, config['max_num_epochs']):

            epoch_start_time = time.time()
            train(epoch)
            epoch_vld_loss, epoch_vld_acc = evaluate(vld_dataloader)
            if nni_search:
                nni.report_intermediate_result({
                    'vld_avg_loss': epoch_vld_loss,
                    'vld_acc': epoch_vld_acc,
                    'default': epoch_vld_acc,
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
                best_model = copy.deepcopy(model)
            elif epoch - best_epoch >= config['early_stopping_patience']:
                print('exiting from training early for early stopping ... ')
                break

            if math.isnan(epoch_vld_loss):
                print('validation loss gets to NaN; exiting current '
                      'invalid trail ...')
                if nni_search:
                    nni.report_final_result({'default': 0})
                sys.exit()
        break

    except RuntimeError as e:
        # CUDA memory or other errors from the training/evaluation process
        traceback.print_exc()
        if nni_search:
            nni.report_final_result({'default': 0})
        sys.exit()

    except KeyboardInterrupt:
        print('exiting from training early for KeyboardInterrupt ... ')
        break

# evaluate the model on the test set
if best_model:

    model = best_model
    tst_start_time = time.time()
    tst_loss, tst_acc = evaluate(tst_dataloader, test=True)
    tst_time_in_sec = time.time() - tst_start_time

    if nni_search:
        nni.report_final_result({
            'tst_loss': tst_loss,
            'tst_acc': tst_acc,
            'default': tst_acc,
        })

    print('=' * 80)
    print(
        f'| end of training '
        f'| test time {tst_time_in_sec:>5.0f} s '
        f'| test loss {tst_loss:5.4f} '
        f'| test accuracy {(tst_acc * 100):3.2f}% '
    )
    print('=' * 80)
