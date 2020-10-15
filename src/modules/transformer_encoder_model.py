"""
File Name:          transformer.py
Project:            bioseq-learning

File Description:

"""
from typing import Optional
from collections import OrderedDict

import torch
import torch.nn as nn

from src.modules.positional_encoding import PositionalEncoding


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
            xfmr_enc_layer_norm: bool,
            xfmr_enc_num_layers: int,
    ):
        super(TransformerEncoderModel, self).__init__()

        self._seq_len: int = seq_len
        self._pos_enc: bool = pos_enc
        self._xfmr_enc_layer_norm: bool = xfmr_enc_layer_norm
        self._xfmr_enc_num_layers: int = xfmr_enc_num_layers

        # TODO: add k-mer embedding option
        _layers = OrderedDict()
        _layers['emb'] = nn.Embedding(
            num_embeddings=num_tokens,
            embedding_dim=emb_dim,
            padding_idx=padding_index,
        )
        if pos_enc:
            _layers['pos_enc'] = PositionalEncoding(
                seq_len=seq_len,
                emb_dim=emb_dim,
                dropout=pos_enc_dropout,
                emb_scale=pos_enc_emb_scale,
            )
        for _i in range(1, xfmr_enc_num_layers + 1):
            _layers[f'xfmr_enc_layer_{_i}'] = nn.TransformerEncoderLayer(
                d_model=emb_dim,
                nhead=xfmr_enc_layer_num_attn_heads,
                dim_feedforward=xfmr_enc_layer_feedforward_dim,
                activation=xfmr_enc_layer_activation,
                dropout=xfmr_enc_layer_dropout,
            )
        if xfmr_enc_layer_norm:
            _layers['xfmr_enc_layer_norm'] = nn.LayerNorm(
                normalized_shape=[seq_len, batch_size, emb_dim]
            )
        _layers['dec'] = nn.Linear(emb_dim, num_tokens)
        self._layers = nn.Sequential(_layers)

        # self._init_weights()

    def _init_weights(self):
        # the initialization methods and parameters are purely random
        nn.init.normal_(getattr(self._layers, 'emb').weight, mean=0.0, std=1.0)
        nn.init.normal_(getattr(self._layers, 'dec').weight, mean=0.0, std=1.0)

    def forward(
            self,
            src: torch.LongTensor,
            attn_mask: Optional[torch.FloatTensor],
            padding_mask: Optional[torch.BoolTensor],
    ) -> torch.FloatTensor:

        _tmp = getattr(self._layers, 'emb')(src)
        if self._pos_enc:
            _tmp = getattr(self._layers, 'pos_enc')(_tmp)
        for _i in range(1, self._xfmr_enc_num_layers + 1):
            _tmp = getattr(self._layers, f'xfmr_enc_layer_{_i}')(
                src=_tmp,
                src_mask=attn_mask,
                src_key_padding_mask=padding_mask,
            )
        return getattr(self._layers, 'dec')(_tmp)
