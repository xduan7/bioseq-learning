"""
File Name:          transformer.py
Project:            bioseq-learning

File Description:

"""
import logging
import random
from collections import OrderedDict
from typing import Iterable, Optional, Tuple

import torch
import torch.nn as nn

from src.modules.positional_encoding import PositionalEncoding


_LOGGER = logging.getLogger(__name__)

MaskerInput = Tuple[
    torch.FloatTensor,
    Optional[torch.BoolTensor]
]
TransformerInput = Tuple[
    torch.FloatTensor,
    Optional[torch.BoolTensor],
    Optional[torch.BoolTensor],
]


class Masker(nn.Module):
    def __init__(
            self,
            seq_len: int,
            num_masks: int,
            mask_index: int,
            attn_mask: bool,
    ):
        super(Masker, self).__init__()

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
            indexed_seqs: torch.FloatTensor,
    ) -> Tuple[torch.FloatTensor, Optional[torch.BoolTensor]]:

        _masked_indexed_seqs = indexed_seqs.clone()
        _masked_indexed_seqs[:, self.src_mask] = self._mask_index
        _attn_mask = self.attn_mask.data
        return _masked_indexed_seqs, _attn_mask


class Seq2Kmer(nn.Module):
    def __init__(
            self,
            k: int,
            seq_len: int,
            num_tokens: int,
            attn_mask: bool,
            padding_mask: bool,
    ):
        super(Seq2Kmer, self).__init__()

        self._k: int = k
        self._seq_len: int = seq_len
        self._num_tokens: int = num_tokens
        self._attn_mask: bool = attn_mask
        self._padding_mask: bool = padding_mask
        self._kmer_seq_len: int = self._seq_len - self._k + 1

        # if padding mask is still requested for k-mer sequence, which is
        # actually not well-defined, return an all-False tensor indicating
        # that attention will be applied everywhere on the k-mer sequence
        if (self._k > 1) and self._attn_mask:
            _warning_msg = \
                f'Attention mask on k-mer (k={self._k}) genome sequence ' \
                f'with masked nucleotide bases is not well-defined. Will ' \
                f'pass-in all \'False\' tensor indicating that the ' \
                f'model shall attend all positions of the sequence.'
            _LOGGER.warning(_warning_msg)
            self.attn_mask: nn.Parameter = nn.Parameter(
                torch.zeros(
                    size=(self._kmer_seq_len, self._kmer_seq_len),
                    dtype=torch.bool,
                ),
                requires_grad=False,
            )

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:

        if self._k == 1:
            return xfmr_input

        indexed_seq, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]

        _indexed_kmer_seq = torch.zeros(
            (indexed_seq.shape[0], self._kmer_seq_len),
            dtype=indexed_seq.dtype,
        ).to(indexed_seq.device)
        for __k in range(self._k):
            _indexed_kmer_seq += \
                indexed_seq[:, __k: self._kmer_seq_len + __k] * \
                (self._num_tokens ** (self._k - __k - 1))

        if self._attn_mask and (attn_mask is not None):
            attn_mask = self.attn_mask.data
        if self._padding_mask and (padding_mask is not None):
            padding_mask = padding_mask[:, :self._kmer_seq_len]

        return _indexed_kmer_seq, attn_mask, padding_mask


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
        # nn.init.normal_(self._emb.weight, mean=0.0, std=1.0)
        for __p in self._emb.parameters():
            if __p.dim() > 1:
                nn.init.orthogonal_(__p)

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
            pos_enc_trainable: bool,
    ):
        super(_PositionalEncoding, self).__init__()
        self._pos_enc = PositionalEncoding(
            seq_len=seq_len,
            emb_dim=emb_dim,
            dropout=pos_enc_dropout,
            emb_scale=pos_enc_emb_scale,
            trainable=pos_enc_trainable,
        ) if pos_enc else None

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]
        _tmp = src.transpose(0, 1).contiguous()
        _tmp = self._pos_enc(_tmp) if self._pos_enc else _tmp
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
            xfmr_padding_mask: bool,
    ):
        super(_TransformerEncoderLayer, self).__init__()
        self._attn_mask: bool = xfmr_attn_mask
        self._padding_mask: bool = xfmr_padding_mask
        self._xfmr_enc_layer = nn.TransformerEncoderLayer(
            d_model=emb_dim,
            nhead=xfmr_enc_layer_num_attn_heads,
            dim_feedforward=xfmr_enc_layer_feedforward_dim,
            activation=xfmr_enc_layer_activation,
            dropout=xfmr_enc_layer_dropout,
        )
        self._init_weights()

    def _init_weights(self):
        for __p in self._xfmr_enc_layer.parameters():
            if __p.dim() > 1:
                nn.init.orthogonal_(__p)
                # nn.init.normal_(__p, mean=0.0, std=1.0)

    def forward(self, xfmr_input: TransformerInput) -> TransformerInput:
        src, attn_mask, padding_mask = \
            xfmr_input[0], xfmr_input[1], xfmr_input[2]

        _tmp = self._xfmr_enc_layer(
            src=src,
            src_mask=attn_mask if self._attn_mask else None,
            src_key_padding_mask=padding_mask if self._padding_mask else None,
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
        self._dec = nn.Sequential(
            # nn.Linear(emb_dim, emb_dim),
            nn.Linear(emb_dim, num_tokens),
        )
        self._init_weights()

    def _init_weights(self):
        for __p in self._dec.parameters():
            if __p.dim() > 1:
                nn.init.orthogonal_(__p)
                # nn.init.normal_(__p, mean=0.0, std=1.0)

    def forward(self, xfmr_input: TransformerInput) -> torch.Tensor:
        src = xfmr_input[0]
        return self._dec(src).transpose(0, 1).contiguous()


class Kmer2Seq(nn.Module):
    def __init__(
            self,
            k: int,
            seq_len: int,
            num_tokens: int,

    ):
        super(Kmer2Seq, self).__init__()

        self._k: int = k
        self._seq_len: int = seq_len
        self._num_tokens: int = num_tokens
        self._kmer_seq_len: int = self._seq_len - self._k + 1

    def forward(self, indexed_kmer_seq: torch.Tensor) -> torch.Tensor:

        if self._k == 1:
            return indexed_kmer_seq

        _indexed_kmer_seq = indexed_kmer_seq.clone()
        _indexed_seq = torch.zeros(
            (_indexed_kmer_seq.shape[0], self._seq_len),
            dtype=_indexed_kmer_seq.dtype,
        ).to(indexed_kmer_seq.device)
        for __k in range(self._k):
            __factor: int = (self._num_tokens ** (self._k - __k - 1))
            if __k == 0:
                _indexed_seq[:, 0: self._kmer_seq_len] += \
                    (_indexed_kmer_seq // __factor)
            else:
                _indexed_seq[:, self._kmer_seq_len + __k - 1] = \
                    (_indexed_kmer_seq[:, -1] // __factor)
            _indexed_kmer_seq = _indexed_kmer_seq % __factor

        return _indexed_seq


def get_transformer_encoder_model(
        num_tokens: int,
        padding_index: int,
        seq_len: int,
        emb_dim: int,
        pos_enc: bool,
        pos_enc_dropout: float,
        pos_enc_trainable: bool,
        pos_enc_emb_scale: float,
        xfmr_enc_layer_num_attn_heads: int,
        xfmr_enc_layer_feedforward_dim: int,
        xfmr_enc_layer_activation: str,
        xfmr_enc_layer_dropout: float,
        xfmr_enc_layer_norm: bool,
        xfmr_enc_num_layers: int,
        xfmr_attn_mask: bool,
        xfmr_padding_mask: bool,
) -> nn.Sequential:
    # this function returns a model suitable for GPipe parallelization
    # which means that the input of the model should be a single tensor or a
    # tuple of them, with batch size as the first dimension
    # otherwise, if this condition cannot be satisfied, then GPipe should be
    # used with num_chunks=1, which means virtually no pipeline and GPipe is
    # only good for load-balancing between multiple GPUs
    _layers = OrderedDict()
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
        pos_enc_trainable=pos_enc_trainable,
    )
    for _i in range(1, xfmr_enc_num_layers + 1):
        _layers[f'xfmr_enc_layer_{_i}'] = _TransformerEncoderLayer(
            emb_dim=emb_dim,
            xfmr_enc_layer_num_attn_heads=xfmr_enc_layer_num_attn_heads,
            xfmr_enc_layer_feedforward_dim=xfmr_enc_layer_feedforward_dim,
            xfmr_enc_layer_activation=xfmr_enc_layer_activation,
            xfmr_enc_layer_dropout=xfmr_enc_layer_dropout,
            xfmr_attn_mask=xfmr_attn_mask,
            xfmr_padding_mask=xfmr_padding_mask,
        )
    _layers['xfmr_enc_layer_norm'] = _LayerNorm(
        emb_dim=emb_dim,
        xfmr_enc_layer_norm=xfmr_enc_layer_norm,
    )
    _layers['dec'] = _Decoder(num_tokens=num_tokens, emb_dim=emb_dim)
    return nn.Sequential(_layers)


def visualize_embedding(
        model: nn.Sequential,
        k: int,
        num_tokens: int,
):
    """TODO: move it somewhere else

    :param model:
    :type model:
    :param k:
    :type k:
    :param num_tokens:
    :type num_tokens:
    :return:
    :rtype:
    """
    import os
    import pandas as pd
    import seaborn as sns
    from copy import deepcopy
    from itertools import product
    import matplotlib.pyplot as plt

    from umap import UMAP
    from Bio.Seq import Seq

    from datasets.genome_dataset import NUCLEOTIDE_CHAR_INDEX_DICT

    _model = deepcopy(model).to('cpu')
    _seq2kmer = Seq2Kmer(
        k=k,
        seq_len=3,
        num_tokens=num_tokens,
        attn_mask=False,
        padding_mask=False,
    )
    _emb = _model.emb

    _codons = [''.join(__p) for __p in product('atgc', repeat=3)]
    _codon_names = [f'{Seq(__c).translate()}({__c})' for __c in _codons]
    _codon_index_seqs = [
        [NUCLEOTIDE_CHAR_INDEX_DICT[__b] for __b in __c] for __c in _codons
    ]

    _codon_index_seqs = torch.LongTensor(_codon_index_seqs).to()
    _codon_kmer_seq, _, _ = _seq2kmer((_codon_index_seqs, None, None))
    _emb_kmer_seq, _, _ = _emb((_codon_kmer_seq, None, None))
    _emb_kmer_seq = _emb_kmer_seq.detach().view(len(_codons), -1).numpy()


    fig_dir_path = './docs/images/embedding/trainable'
    os.makedirs(fig_dir_path, exist_ok=True)

    for __n_neighbors in [2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32]:
        for __min_dist in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
            umap = UMAP(
                n_neighbors=__n_neighbors,
                min_dist=__min_dist,
                init='random',
            )
            _emb_kmer_coord = umap.fit_transform(_emb_kmer_seq)
            _emb_kmer_coord_df = pd.DataFrame(
                _emb_kmer_coord,
                index=_codon_names,
                columns=['x', 'y'],
            )
            _emb_kmer_coord_df['animo acid'] = \
                _emb_kmer_coord_df.index.str.split('(').str.get(0)

            plt.figure(figsize=(8, 6))
            sns.scatterplot(
                data=_emb_kmer_coord_df,
                x='x', y='y',
                hue='animo acid',
            )
            plt.legend(bbox_to_anchor=(1, 1), loc='upper left')

            fig_name = \
                f'umap_n_neighbors_{__n_neighbors:06d}_' \
                f'umap_min_dist_{__min_dist:.3f}.png'
            plt.savefig(os.path.join(fig_dir_path, fig_name))
