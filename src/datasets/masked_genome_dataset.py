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


MASK_CHAR: str = '*'
PADDING_CHAR: str = '-'
NUCLEOTIDE_CHAR_SET: Set[str] = {'a', 't', 'g', 'c'}
NUCLEOTIDE_CHAR_INDEX_DICT: Dict[str, int] = {
    MASK_CHAR: 5,
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

        # dict that maps (genome id + contig id) to contig sequence
        self._genome_contig_seq_dict: Dict[str, str] = {}

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
                # note that this does not work for really short contigs
                # e.g. when seq_len > len(contig). in this case,
                # the number of paddings will exceed the given maximum
                _contig_seq_len: int = len(_contig_seq_rec.seq)
                if _contig_seq_len < self._seq_len:
                    _warning_msg = \
                        f'The length of contig {_contig_seq_rec.id} from ' \
                        f'genome {_genome_id} is {_contig_seq_len}, ' \
                        f'smaller compared to the intended dataset ' \
                        f'sequence length {self._seq_len}. ' \
                        f'Ignoring the contig ...'
                    _LOGGER.warning(_warning_msg)
                self._genome_contig_seq_dict[_genome_contig_id] = \
                    max_num_paddings * PADDING_CHAR + \
                    str(_contig_seq_rec.seq).lower() + \
                    max_num_paddings * PADDING_CHAR

        # dict that maps dataset index (accessible from dataloader) to
        # genome id + contig id, and the starting position of the sequence
        # check __getitem__ function for the usage details
        _index_counter = 0
        self._index_genome_contig_id_pos_dict: Dict[int, Tuple[str, int]] = {}
        for _genome_contig_id, _padded_contig_seq in \
                self._genome_contig_seq_dict.items():
            for _pos in range(len(_padded_contig_seq) - self._seq_len + 1):
                self._index_genome_contig_id_pos_dict[_index_counter] = \
                    (_genome_contig_id, _pos)
                _index_counter += 1

        self.__len = len(self._index_genome_contig_id_pos_dict)

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
        _genome_contig_id, _pos = \
            self._index_genome_contig_id_pos_dict[index]
        _genome_contig_seq: str = \
            self._genome_contig_seq_dict[_genome_contig_id]
        _seq: str = _genome_contig_seq[_pos: _pos + self._seq_len]
        _seq_char_list: List[str] = list(_seq)

        # mask the nucleotide randomly
        _nucleotide_indices = [
            _i for _i, _c in enumerate(_seq_char_list)
            if _c in NUCLEOTIDE_CHAR_SET]
        _masked_indices: List[int] = \
            random.sample(_nucleotide_indices, self._num_masks)

        # # masked the (not indexed) sequence
        # for _i in _masked_indices:
        #     _seq_char_list[_i] = MASK_CHAR
        # _masked_seq = ''.join(_seq_char_list)

        _indexed_seq = torch.LongTensor(
            list(map(NUCLEOTIDE_CHAR_INDEX_DICT.get, _seq)))

        _mask = torch.zeros_like(_indexed_seq, dtype=torch.float)
        _mask.scatter_(0, torch.LongTensor(_masked_indices), float('-inf'))

        return _indexed_seq, _mask
