"""
File Name:          sequence_mask.py
Project:            bioseq-learning

File Description:

"""
import random
from typing import Tuple, Iterable, Optional

import torch


class SequenceMask():

    def __init__(
            self,
            seq_len: int,
            num_masks: int,
            init_masked_indices: Optional[Iterable[int]] = None,
    ):
        """
        TODO: docstring
        :param seq_len:
        :type seq_len:
        :param init_masked_indices:
        :type init_masked_indices:
        :return:
        :rtype:
        """
        assert 0 < num_masks < seq_len

        self._seq_len: int = seq_len
        self._num_masks: int = num_masks

        # generate the default masks for both attention and source sequences
        self._attn_mask: torch.FloatTensor
        self._src_mask: torch.BoolTensor
        self._attn_mask, self._src_mask = self.update(init_masked_indices)

    def update(
            self,
            masked_indices: Optional[Iterable[int]] = None,
    ) -> Tuple[torch.FloatTensor, torch.BoolTensor]:
        """
        TODO: docstring
        :param masked_indices:
        :type masked_indices:
        :return:
        :rtype:
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

    def get_attn_mask(self):
        """
        TODO: docstring
        :return:
        :rtype:
        """
        return self._attn_mask

    def get_src_mask(self):
        """
        TODO: docstring
        :return:
        :rtype:
        """
        return self._src_mask
