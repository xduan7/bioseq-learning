"""
File Name:          positional_encoding.py
Project:            bioseq-learning

File Description:

    Official PyTorch implementation of positional encoding, copied from
    the transformer tutorial with some changes on the comments:
    https://pytorch.org/tutorials/beginner/transformer_tutorial.html

    My own previous implementation which is functionally the same:
    https://github.com/xduan7/DLTM/blob/master/networks/transformer/positional_encoder.py

    The idea of positional encoding is from the paper [Attention Is All You
    Need](https://arxiv.org/abs/1706.03762)

"""
import math
import torch
import torch.nn as nn


class PositionalEncoding(nn.Module):
    """positional encoding module proposed in "Attention Is All You Need"
    """

    def __init__(
            self,
            seq_len: int,
            emb_dim: int,
            dropout: float = 0.1,
            emb_scale: float = 1.0,
            trainable: bool = False,
    ):
        """constructor for positional encoding

        Create a positional encoding matrix that contains the information
        about the relative or absolute position of the tokens in the
        sequence. The positional encodings have the same dimension as the
        embeddings, so that the two can be summed. Here, we use sine and
        cosine functions of different frequencies.
        .. math::
            \text{PosEnc}(pos, 2i) = sin(pos/10000^(2i/emb_dim))
            \text{PosEnc}(pos, 2i+1) = cos(pos/10000^(2i/emb_dim))
            \text{where pos is the word position and i is the embed idx)

        :param seq_len: (maximum) length of the input sequence
        :type seq_len: int
        :param emb_dim: embedding dimension for sequence model; must be an
        even number for both sin and cos positional terms
        :type emb_dim: int
        :param dropout: dropout rate (after adding positional encoding); de
        :type dropout: float
        :param emb_scale: scaling factor for embedded vectors to prevent
        them from getting diminishing after adding positional encodings
        :type emb_scale: float
        :param trainable: TODO
        :type trainable: bool
        :return: None
        """
        super(PositionalEncoding, self).__init__()

        # embedding dimension must be even for both sin and cos terms
        assert emb_dim % 2 == 0
        self._seq_len: int = seq_len
        self._emb_dim: int = emb_dim

        # scaling embedded vector to prevent it from getting diminished after
        # adding positional encoding
        self._emb_scale = emb_scale
        self._dropout = nn.Dropout(p=dropout)

        positional_encoding = torch.zeros(seq_len, emb_dim)
        _positions = torch.arange(0, seq_len, dtype=torch.float)
        _positions = _positions.unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, emb_dim, 2).float()
            * (-math.log(10000.0) / emb_dim)
        )
        positional_encoding[:, 0::2] = torch.sin(_positions * div_term)
        positional_encoding[:, 1::2] = torch.cos(_positions * div_term)
        positional_encoding = \
            positional_encoding.unsqueeze(0).transpose(0, 1)

        if trainable:
            self.register_parameter(
                'positional_encoding',
                nn.Parameter(positional_encoding, requires_grad=True),
            )
        else:
            # use register_buffer (returns tensor) instead of nn.Parameter()
            # so that the matrix is not learnable/trainable
            self.register_buffer('positional_encoding', positional_encoding)

    def forward(
            self,
            x: torch.Tensor,
    ) -> torch.Tensor:
        """add position encoding information to scaled embedded vector of a
        batch of sequences, and return after applying the dropout.

        :param x: input embedded sequence tensor
        :type x: torch.Tensor in the shape of (seq_len, batch_size, emb_dim)
        :return: torch.Tensor in the shape of (seq_len, batch_size, emb_dim)
        """
        x = x * self._emb_scale + self.positional_encoding[:x.size(0), :]
        return self._dropout(x)

    def plot(self, file_path: str):
        """TODO

        :param file_path:
        :type file_path:
        :return:
        :rtype:
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from torch.autograd import Variable

        plt.figure(figsize=(48, 8))
        _device = self.positional_encoding.device
        _input = Variable(torch.zeros(self._seq_len, 1, self._emb_dim))
        _output = self.forward(_input.to(_device)).cpu()
        plt.plot(np.arange(self._seq_len), _output[:, 0, 0::64].data.numpy())

        plt.savefig(file_path)
        plt.close('all')
