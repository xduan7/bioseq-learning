"""
File Name:          masked_genome_dataset.py
Project:            bioseq-learning-cd

File Description:

"""
import os
import random
import logging
from typing import Tuple, List, Set, Dict, Iterable, Union

import torch
import numpy as np
from Bio import SeqIO, SeqRecord
from torch.utils.data import Dataset


# MASK_CHAR: str = '*'
PADDING_CHAR: str = '-'
NUCLEOTIDE_CHAR_SET: Set[str] = {'a', 't', 'g', 'c'}
NUCLEOTIDE_CHAR_INDEX_DICT: Dict[str, int] = {
    # MASK_CHAR: 5,
    PADDING_CHAR: 0,
    'a': 1,
    't': 2,
    'g': 3,
    'c': 4,
}
NUCLEOTIDE_CHAR_VOCAB_SIZE: int = len(NUCLEOTIDE_CHAR_INDEX_DICT)

_LOGGER = logging.getLogger(__name__)


class MaskedGenomeDataset(Dataset):
    """basic genome dataset for (dynamic) masked language model training
    """
    def __init__(
            self,
            genome_dir_paths: Iterable[str],
            seq_len: int,
            num_masks: Union[int, float],
            max_num_paddings: int,

    ):
        """
        :param genome_dir_paths: an iterable of paths to individual genomes
        that are processed (check process_genome.py for details)
        :type genome_dir_paths: iterable of strings
        :param seq_len: the length of genome sequences of the dataset after
        segmentation and padding
        :type seq_len: int
        :param num_masks: number or percentage of masks in applied in a
        genome sequence; note that paddings are not applicable for masking
        :type num_masks: int or float
        :param max_num_paddings: maximum number of padding characters in a
        genome sequence
        :type max_num_paddings: int
        """

        # get the actual number of masks if the argument is given as a float
        num_masks: int = int(np.round(num_masks * seq_len)) \
            if isinstance(num_masks, float) else num_masks

        # sanity check for the sequence length, paddings, and masks
        # TODO: add ValueError or some other exceptions
        assert 0 < max_num_paddings < seq_len
        assert 0 < num_masks < (seq_len - max_num_paddings)

        self._seq_len: int = seq_len
        self._num_masks: int = num_masks
        self._max_num_paddings: int = max_num_paddings
        self.__len = 0

        # dict that maps (genome id + contig id) to contig sequence
        self._genome_contig_seq_dict: Dict[str, str] = {}
        # data structure that stores a iterable of tuples that contains
        # (genome id + contig id, number of samples/segmented sequences)
        _genome_contig_num_seqs_list: List[Tuple[str, int]] = []
        for _genome_dir_path in genome_dir_paths:
            _genome_id: str = os.path.basename(_genome_dir_path.rstrip('/'))
            _genome_contig_seq_dir_path: str = os.path.join(
                _genome_dir_path, 'contigs')
            # _genome_feature_dir_path: str = os.path.join(
            #     _genome_dir_path, 'features')
            for _contig_seq_file_name in \
                    os.listdir(_genome_contig_seq_dir_path):
                _contig_seq_path = os.path.join(
                    _genome_contig_seq_dir_path,
                    _contig_seq_file_name,
                )
                with open(_contig_seq_path, 'r') as _fh:
                    _contig_seq_rec: SeqRecord = \
                        next(SeqIO.parse(_fh, 'fasta'))
                _genome_contig_id: str = \
                    f'{_genome_id}/{_contig_seq_rec.id}'
                _padded_seq: str = \
                    max_num_paddings * PADDING_CHAR + \
                    str(_contig_seq_rec.seq).lower() + \
                    max_num_paddings * PADDING_CHAR
                # note that this does not work for really short contigs
                # e.g. when 'seq_len + paddings > len(contig)'. in this case,
                # the number of paddings will exceed the given maximum
                if len(_padded_seq) < self._seq_len:
                    _warning_msg = \
                        f'The length of contig {_contig_seq_rec.id} from ' \
                        f'genome {_genome_id} is {len(_padded_seq)} with ' \
                        f'paddings on both ends, smaller compared to the ' \
                        f'intended dataset sequence length {self._seq_len}' \
                        f'. Ignoring the contig ...'
                    _LOGGER.warning(_warning_msg)
                self._genome_contig_seq_dict[_genome_contig_id] = _padded_seq
                _num_seqs: int = len(_padded_seq) - self._seq_len + 1
                _genome_contig_num_seqs_list.append(
                    (_genome_contig_id, _num_seqs))
                self.__len += _num_seqs

        # convert list to tuple for faster access
        self._genome_contig_num_seqs_tuple: Tuple[Tuple[str, int], ...] = \
            tuple(_genome_contig_num_seqs_list)


    def __len__(self) -> int:
        return self.__len

    def __getitem__(
            self,
            index: int,
    ):
        """
        :param index: index to the genome sequence in [0, len(self)]
        :type index: int
        :return: a tuple that contains a indexed genome sequence in
        LongTensor amd a mask tensor of the same size, but in FloatTensor
        where float('-inf') represents a mask
        :rtype: tuple of tensors
        """
        if index < 0:
            index = index + self.__len

        _seq, _seq_char_list = '', []
        for _genome_contig_id, _num_seqs in self._genome_contig_num_seqs_tuple:
            if index >= _num_seqs:
                index = index - _num_seqs
            else:
                _genome_contig_seq: str = \
                    self._genome_contig_seq_dict[_genome_contig_id]
                _seq: str = _genome_contig_seq[index: index + self._seq_len]
                _seq_char_list: List[str] = list(_seq)
                break

        # convert the nucleotide sequence to the indexed (numeric) sequence
        _indexed_seq = torch.LongTensor(
            list(map(NUCLEOTIDE_CHAR_INDEX_DICT.get, _seq)))

        # mask the nucleotide at random location (not on paddings though)
        _nucleotide_indices = [
            _i for _i, _c in enumerate(_seq_char_list)
            if _c in NUCLEOTIDE_CHAR_SET]
        _masked_indices: List[int] = \
            random.sample(_nucleotide_indices, self._num_masks)
        _mask = torch.zeros_like(_indexed_seq, dtype=torch.float)
        _mask.scatter_(0, torch.LongTensor(_masked_indices), float('-inf'))

        return _indexed_seq, _mask
