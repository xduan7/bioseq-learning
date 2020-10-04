"""
File Name:          transformer.py
Project:            bioseq-learning

File Description:

"""
import copy
import random
from typing import Tuple

import torch
import torch.nn as nn

from src.modules import PositionalEncoding


class TransformerEncoderModel(nn.Module):
    """container module for a complete end-to-end transformer encoder model,
    with a embedding layer, (possibly) a position encoding layer,
    and a linear decoding layer at the end
    """
    def __init__(
            self,
            num_tokens: int,
            padding_index: int,
            seq_len: int,
            batch_size: int,
            emb_dim: int,
            pos_enc: bool,
            pos_enc_dropout: float,
            pos_enc_emb_scale: float,
            xfmr_enc_layer_num_attn_heads: int,
            xfmr_enc_layer_feedforward_dim: int,
            xfmr_enc_layer_activation: str,
            xfmr_enc_layer_dropout: float,
            xfmr_enc_num_layers: int,
            xfmr_enc_norm: bool,
    ):
        super(TransformerEncoderModel, self).__init__()

        self._seq_len = seq_len

        # TODO: add k-mer embedding option
        self._emb = nn.Embedding(
            num_embeddings=num_tokens,
            embedding_dim=emb_dim,
            padding_idx=padding_index,
        )
        self._pos_enc = PositionalEncoding(
            seq_len=seq_len,
            emb_dim=emb_dim,
            dropout=pos_enc_dropout,
            emb_scale=pos_enc_emb_scale,
        ) if pos_enc else None
        _xfmr_enc_layer = nn.TransformerEncoderLayer(
            d_model=emb_dim,
            nhead=xfmr_enc_layer_num_attn_heads,
            dim_feedforward=xfmr_enc_layer_feedforward_dim,
            activation=xfmr_enc_layer_activation,
            dropout=xfmr_enc_layer_dropout,
        )
        _xfmr_enc_norm = nn.LayerNorm(
            normalized_shape=[seq_len, batch_size, emb_dim]
        ) if xfmr_enc_norm else None
        self._xfmr_enc = nn.TransformerEncoder(
            encoder_layer=_xfmr_enc_layer,
            num_layers=xfmr_enc_num_layers,
            norm=_xfmr_enc_norm,
        )
        self._dec = nn.Linear(emb_dim, num_tokens)

        self._init_weights()

    def _init_weights(self):
        # the initialization methods and parameters are purely random
        nn.init.normal_(self._emb.weight, mean=0.0, std=1.0)
        nn.init.normal_(self._dec.weight, mean=0.0, std=1.0)

    def forward(
            self,
            src: torch.LongTensor,
            attn_mask: torch.FloatTensor,
            padding_mask: torch.BoolTensor,
    ) -> torch.FloatTensor:

        # dynamically generate a mask of size (seq_len, seq_len)
        # indicates the accessibility of the tokens for the ith sample
        # _masked_indices = random.sample(range(self._seq_len), num_masks)
        # _mask = torch.zeros(size=(self._seq_len, ), dtype=torch.float)
        # _mask.scatter_(0, torch.LongTensor(_masked_indices), float('-inf'))
        # _mask = _mask.repeat(self._seq_len).view(-1, self._seq_len)
        # _mask = _mask.to(src.device)
        #
        # _tmp = copy.deepcopy(src)
        # _tmp[_mask[0].bool()] = 0

        _tmp = self._emb(src)
        if self._pos_enc:
            _tmp = self._pos_enc(_tmp)
        _tmp = self._xfmr_enc(
            src=_tmp,
            mask=attn_mask,
            src_key_padding_mask=padding_mask,
        )
        return self._dec(_tmp)
