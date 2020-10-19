"""
File Name:          profile_sequential_module.py
Project:            bioseq-learning

File Description:

"""
from typing import List, Sequence, Union

import torch
import torch.nn as nn
from torch import Tensor
from torchgpipe.microbatch import Batch
from torchgpipe.balance.profile import layerwise_sandbox, detach


def profile_sequential_module(
        module: nn.Sequential,
        input: Union[Tensor, Sequence[Tensor]],
        chunks: int,
        param_scale: float,
        device: torch.device,
) -> List[int]:
    """similar to 'profile_sizes' function in torchgpipe, but instead of
    passing in a batch of size 1, it passes in a whole batch for more
    accurate estimate of the sizes; moreover, it fixed the issue with
    negative memory allocation for latent variables

    reference: torchgpipe.balance.profile.profile_sizes

    :param module: pytorch sequential module to be profiled
    :type module: nn.Sequential
    :param input: input tensor or a sequence (will be cast to tuple) of tensors
    :type input: Union[Tensor, Sequence[Tensor]]
    :param chunks: number of chunks for a single batch specified in GPipe
    :type chunks: int
    :param param_scale: scaling factor for parameters (SGD: 2-3, Adam: 4-5,
    etc.); check GPipe doc for more details
    more details
    :type param_scale: float
    :param device: device for size profiling run; must be GPU
    :type device: torch.device
    :return: list of integers representing the sizes of all the layers in
    sequential model in bytes
    :rtype: List[int]
    """
    if device.type != 'cuda':
        raise ValueError('require CUDA device for size profiler supports '
                         'only CUDA device')

    # cast everything in the batch into a tuple of tensors if the given
    # input is a sequence of tensors
    _batch = Batch(input) if isinstance(input, Tensor) else \
        Batch(tuple([_i.detach().to(device) for _i in input]))
    _layer_sizes_in_byte: List[int] = []

    for layer in layerwise_sandbox(module, device):
        detach(_batch)

        # Detect memory usage at forward.
        _memory_before = torch.cuda.memory_allocated(device)
        _batch = _batch.call(layer)
        _memory_after = torch.cuda.memory_allocated(device)
        _latent_size = max(0, _memory_after - _memory_before)

        # Analyze size of parameters.
        param_size = sum(
            p.storage().size() * p.storage().element_size()
            for p in layer.parameters())

        # Combine size of parameters and activations with normalize
        # scales.
        _size = _latent_size / chunks + param_size * param_scale
        _layer_sizes_in_byte.append(int(_size))

    return _layer_sizes_in_byte
