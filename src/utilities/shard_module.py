"""
File Name:          shard_module.py
Project:            bioseq-learning

File Description:

"""
import logging
from typing import List, Sequence, Union

import torch
import torch.nn as nn
from torch import Tensor, device
from torchgpipe import GPipe
from torchgpipe.balance.blockpartition import solve
from src.utilities.profile_sequential_module import \
    profile_sequential_module


_LOGGER = logging.getLogger(__name__)


def shard_module(
        module: nn.Sequential,
        input_sample: Union[Tensor, Sequence[Tensor]],
        devices: Sequence[device],
        num_chunks: int,
) -> Union[GPipe, nn.Sequential]:

    # if more than 1 devices is necessary for the run ...
    # gpipe model on 1 gpu is much slower compared to ordinary model
    if len(devices) > 1:

        # get the sizes by layer in the sequential model
        # merely an rough estimation with a forward batch without optimizer
        _module_layer_sizes_in_byte: List[int] = profile_sequential_module(
            module=module,
            input=input_sample,
            chunks=1,
            param_scale=5.0,
            device=devices[0],
        )

        # scale the sizes by 2 (forward and backward) to play it safe
        _module_layer_sizes_in_byte = \
            [__m * 2 for __m in _module_layer_sizes_in_byte]

        # get the minimal necessary number of devices
        _num_devices: int = 0
        _module_size_in_byte: int = sum(_module_layer_sizes_in_byte)
        _gpu_sizes_in_byte: List[int] = [
            torch.cuda.get_device_properties(__d).total_memory
            for __d in devices
        ]

        for __n in range(1, len(devices) + 1):
            if sum(_gpu_sizes_in_byte[:__n]) > _module_size_in_byte:
                _num_devices = __n
                break

        if _num_devices == 0:
            _warning_msg = \
                f'PyTorch module with the estimated size of ' \
                f'{_module_size_in_byte / (1024**3):.1f} Gb exceeds the ' \
                f'combined memory capacity of all computation devices ' \
                f'({devices}). Proceeding with caution ...'
            _LOGGER.warning(_warning_msg)
        elif _num_devices < len(devices):
            _warning_msg = \
                f'Using only {_num_devices} device(s) (' \
                f'{devices[:_num_devices]}) out of {len(devices)} ' \
                f'for the current trial.'
            _LOGGER.warning(_warning_msg)
            devices = devices[:_num_devices]
        else:
            # using all the available devices ...
            pass

        # cast the model to gpipe if necessary (too big for a single GPU)
        if len(devices) > 1:
            _balance = [len(__b) for __b in solve(
                _module_layer_sizes_in_byte,
                partitions=len(devices),
            )]
            # TODO: maybe not using gpipe?
            # since we are not actually chunking the batches ...
            module = GPipe(
                module=module,
                balance=_balance,
                devices=devices,
                chunks=num_chunks,  # subject to change
                checkpoint='never',
            )
        else:
            # if the model can fit into a single GPU
            module = module.to(devices[0])
    else:
        # if the number of available devices is 1
        module = module.to(devices[0])

    return module, devices