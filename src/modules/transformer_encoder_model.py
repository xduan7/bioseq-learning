"""
File Name:          transformer.py
Project:            bioseq-learning

File Description:

"""
import random
from collections import OrderedDict
from typing import Iterable, Optional, Tuple

import torch
import torch.nn as nn

from src.modules.positional_encoding import PositionalEncoding


MaskerInput = Tuple[
    torch.FloatTensor,
    Optional[torch.BoolTensor]
]
TransformerInput = Tuple[
    torch.FloatTensor,
    Optional[torch.BoolTensor],
    Optional[torch.BoolTensor],
]


class _Masker(nn.Module):
    def __init__(
            self,
            seq_len: int,
            num_masks: int,
            mask_index: int,
            attn_mask: bool,
    ):
        super(_Masker, self).__init__()

        assert 0 < num_masks < seq_len
        self._seq_len: int = seq_len
        self._num_masks: int = num_masks
        self._mask_index: int = mask_index
        self._attn_mask: bool = attn_mask

        self.attn_mask: nn.Parameter = nn.Parameter(
            torch.zeros(
                size=(self._seq_len, self._seq_len),
                dtype=torch.bool,
            ) if self._attn_mask else torch.BoolTensor([]),
            requires_grad=False,
        )
        self.src_mask: nn.Parameter = nn.Parameter(
            torch.zeros(
                size=(self._seq_len, ),
                dtype=torch.bool,
            ),
            requires_grad=False,
        )
        self.update()

        # tensors for input/target storage, used for accuracy calculation
        # self.curr_masked_indexed_seqs: Optional[torch.LongTensor] = None

    def update(
            self,
            num_masks: Optional[int] = None,
            masked_indices: Optional[Iterable[int]] = None,
    ):
        """update the attention and sequence masks

        :param num_masks: optional number of masks for update
        :type num_masks int or None
        :param masked_indices: optional indices for masking, randomly
        generate a a list of mask indices if it's not given
        :type masked_indices: iterable of integer or None
        :return None
        :rtype None
        """
        if num_masks:
            self._num_masks = num_masks
        __masked_indices = torch.LongTensor(masked_indices) \
            if masked_indices else torch.LongTensor(
            random.sample(range(self._seq_len), self._num_masks))

        # self.attn_mask.data = None
        # self.src_mask.data = None

        __device = self.src_mask.data.device

        if self._attn_mask:
            __attn_mask: torch.BoolTensor = torch.zeros(
                size=(self._seq_len,), dtype=torch.bool, requires_grad=False,
            ).scatter_(
                0, __masked_indices, True
            ).repeat(self._seq_len).view(-1, self._seq_len)
            self.attn_mask.data = __attn_mask.to(__device)

        __src_mask: torch.BoolTensor = torch.zeros(
            size=(self._seq_len,), dtype=torch.bool, requires_grad=False,
        ).scatter_(0, __masked_indices, True)
        self.src_mask.data = __src_mask.to(__device)

    def forward(
            self,
            input: MaskerInput,
    ) -> TransformerInput:

        _indexed_seqs, _padding_mask = input[0], input[1]
        _masked_indexed_seqs = _indexed_seqs.clone()
        _masked_indexed_seqs[:, self.src_mask] = self._mask_index

        # self.curr_masked_indexed_seqs = _masked_indexed_seqs

        _tmp = _masked_indexed_seqs.transpose(0, 1).contiguous()
        _attn_mask = self.attn_mask.data

        return _tmp, _attn_mask, _padding_mask


class _Embedding(nn.Module):
    def __init__(
            self,
            num_tokens: int,
            emb_dim: int,
            padding_index: int,
    ):
        super(_Embedding, self).__init__()
        self._emb = nn.Embedding(
            num_embeddings=num_tokens,
            embedding_dim=emb_dim,
            padding_idx=padding_index,
        )
        self._init_weights()

    def _init_weights(self):
        nn.init.normal_(self._emb.weight, mean=0.0, std=1.0)

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]
        _tmp = self._emb(src)
        return _tmp, attn_mask, padding_mask


class _PositionalEncoding(nn.Module):
    def __init__(
            self,
            seq_len: int,
            emb_dim: int,
            pos_enc: bool,
            pos_enc_dropout: float,
            pos_enc_emb_scale: float,
    ):
        super(_PositionalEncoding, self).__init__()
        self._pos_enc = PositionalEncoding(
            seq_len=seq_len,
            emb_dim=emb_dim,
            dropout=pos_enc_dropout,
            emb_scale=pos_enc_emb_scale,
        ) if pos_enc else None

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]
        _tmp = self._pos_enc(src) if self._pos_enc else src
        return _tmp, attn_mask, padding_mask


class _TransformerEncoderLayer(nn.Module):
    def __init__(
            self,
            emb_dim: int,
            xfmr_enc_layer_num_attn_heads: int,
            xfmr_enc_layer_feedforward_dim: int,
            xfmr_enc_layer_activation: str,
            xfmr_enc_layer_dropout: float,
            xfmr_attn_mask: bool,
    ):
        super(_TransformerEncoderLayer, self).__init__()
        self._attn_mask: bool = xfmr_attn_mask
        self._xfmr_enc_layer = nn.TransformerEncoderLayer(
            d_model=emb_dim,
            nhead=xfmr_enc_layer_num_attn_heads,
            dim_feedforward=xfmr_enc_layer_feedforward_dim,
            activation=xfmr_enc_layer_activation,
            dropout=xfmr_enc_layer_dropout,
        )

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]
        _tmp = self._xfmr_enc_layer(
            src=src,
            src_mask=attn_mask if self._attn_mask else None,
            src_key_padding_mask=padding_mask,
        )
        return _tmp, attn_mask, padding_mask


class _LayerNorm(nn.Module):
    def __init__(
            self,
            emb_dim: int,
            xfmr_enc_layer_norm: bool,
    ):
        super(_LayerNorm, self).__init__()
        self._layer_norm = nn.LayerNorm(emb_dim) if xfmr_enc_layer_norm else None

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]
        _tmp = self._layer_norm(src) if self._layer_norm else src
        return _tmp, attn_mask, padding_mask


class _Decoder(nn.Module):
    def __init__(
            self,
            num_tokens: int,
            emb_dim: int,
    ):
        super(_Decoder, self).__init__()
        self._dec = nn.Linear(emb_dim, num_tokens)
        self._init_weights()

    def _init_weights(self):
        nn.init.normal_(self._dec.weight, mean=0.0, std=1.0)

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src = xfmr_input[0]
        return self._dec(src).transpose(0, 1).contiguous()


def get_transformer_encoder_model(
        num_tokens: int,
        num_masks: int,
        mask_index: int,
        padding_index: int,
        seq_len: int,
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
        xfmr_attn_mask: bool,
) -> nn.Sequential:

    # TODO: add k-mer embedding option
    _layers = OrderedDict()
    _layers['msk'] = _Masker(
        seq_len=seq_len,
        num_masks=num_masks,
        mask_index=mask_index,
        attn_mask=xfmr_attn_mask,
    )
    _layers['emb'] = _Embedding(
        num_tokens=num_tokens,
        emb_dim=emb_dim,
        padding_index=padding_index,
    )
    _layers['pos_enc'] = _PositionalEncoding(
        seq_len=seq_len,
        emb_dim=emb_dim,
        pos_enc=pos_enc,
        pos_enc_dropout=pos_enc_dropout,
        pos_enc_emb_scale=pos_enc_emb_scale,
    )
    for _i in range(1, xfmr_enc_num_layers + 1):
        _layers[f'xfmr_enc_layer_{_i}'] = _TransformerEncoderLayer(
            emb_dim=emb_dim,
            xfmr_enc_layer_num_attn_heads=xfmr_enc_layer_num_attn_heads,
            xfmr_enc_layer_feedforward_dim=xfmr_enc_layer_feedforward_dim,
            xfmr_enc_layer_activation=xfmr_enc_layer_activation,
            xfmr_enc_layer_dropout=xfmr_enc_layer_dropout,
            xfmr_attn_mask=xfmr_attn_mask,
        )
    _layers['xfmr_enc_layer_norm'] = _LayerNorm(
        emb_dim=emb_dim,
        xfmr_enc_layer_norm=xfmr_enc_layer_norm,
    )
    _layers['dec'] = _Decoder(num_tokens=num_tokens, emb_dim=emb_dim)
    return nn.Sequential(_layers)
