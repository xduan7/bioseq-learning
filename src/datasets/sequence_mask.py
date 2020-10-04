"""
File Name:          sequence_mask.py
Project:            bioseq-learning

File Description:

"""
import random
from typing import Tuple, Iterable, Optional

import torch


class SequenceMask:
    """mask generator class for genome sequences
    """
    def __init__(
            self,
            seq_len: int,
            num_masks: int,
            init_masked_indices: Optional[Iterable[int]] = None,
    ):
        """initialize a mask generator

        :param seq_len:
        :type seq_len:
        :param init_masked_indices:
        :type init_masked_indices:
        """
        assert 0 < num_masks < seq_len

        self._seq_len: int = seq_len
        self._num_masks: int = num_masks

        # generate the default masks for both attention and source sequences
        # attention masks is a float 2D tensor of shape seq_len * seq_len,
        # where -inf represents the attention 'black hole'
        # sequence mask is a boolean 1D tensor of shape seq_len, where True
        # represents masks
        self._attn_mask: torch.FloatTensor
        self._src_mask: torch.BoolTensor
        self._attn_mask, self._src_mask = self.update(init_masked_indices)

    def update(
            self,
            masked_indices: Optional[Iterable[int]] = None,
    ) -> Tuple[torch.FloatTensor, torch.BoolTensor]:
        """update the attention and sequence masks

        :param masked_indices: optional indices for masking, randomly
        generate a a list of mask indices if it's not given
        :type masked_indices: iterable of integer or None
        :return: attention masks (float 2D tensor of shape seq_len * seq_len)
        and the sequence mask (boolean 1D tensor of shape seq_len)
        :rtype: a tuple of float tensor and a boolean tensor
        """

        _masked_indices = masked_indices if masked_indices else \
            random.sample(range(self._seq_len), self._num_masks)
        _masked_indices = torch.LongTensor(_masked_indices)

        _attn_mask: torch.FloatTensor = torch.zeros(
            size=(self._seq_len,),
            dtype=torch.float,
        ).scatter_(0, _masked_indices, float('-inf'))
        self._attn_mask: torch.FloatTensor = \
            _attn_mask.repeat(self._seq_len).view(-1, self._seq_len)

        self._src_mask: torch.BoolTensor = torch.zeros(
            size=(self._seq_len,),
            dtype=torch.bool,
        ).scatter_(0, _masked_indices, True)

        return self._attn_mask, self._src_mask

    def get_attn_mask(self) -> torch.FloatTensor:
        """get the attention mask

        :return: attention masks (float 2D tensor of shape seq_len * seq_len)
        :rtype: float tensor
        """
        return self._attn_mask

    def get_src_mask(self) -> torch.BoolTensor:
        """get the sequence mask

        :return: sequence mask (boolean 1D tensor of shape seq_len)
        :rtype: boolean tensor
        """
        return self._src_mask
