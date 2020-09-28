"""
File Name:          genome_dataset.py
Project:            bioseq-learning-cd

File Description:

"""
import os
import logging
from typing import Tuple, List, Set, Dict, Iterable

import torch
from Bio import SeqIO, SeqRecord
from torch.utils.data import Dataset


PADDING_CHAR: str = '-'
NUCLEOTIDE_CHAR_SET: Set[str] = {'a', 't', 'g', 'c'}
NUCLEOTIDE_CHAR_INDEX_DICT: Dict[str, int] = {
    PADDING_CHAR: 0,
    'a': 1,
    't': 2,
    'g': 3,
    'c': 4,
}
NUCLEOTIDE_CHAR_VOCAB_SIZE: int = len(NUCLEOTIDE_CHAR_INDEX_DICT)
PADDING_INDEX: int = NUCLEOTIDE_CHAR_INDEX_DICT[PADDING_CHAR]

_LOGGER = logging.getLogger(__name__)


class GenomeDataset(Dataset):
    """basic genome dataset for (dynamic) masked language model training
    """
    def __init__(
            self,
            genome_dir_paths: Iterable[str],
            seq_len: int,
            max_num_paddings: int,

    ):
        """
        :param genome_dir_paths: an iterable of paths to individual genomes
        that are processed (check process_genome.py for details)
        :type genome_dir_paths: iterable of strings
        :param seq_len: the length of genome sequences of the dataset after
        segmentation and padding
        :type seq_len: int
        :param max_num_paddings: maximum number of padding characters in a
        genome sequence
        :type max_num_paddings: int
        """

        # sanity check for the sequence length and number of paddings
        assert 0 < max_num_paddings < seq_len

        self._seq_len: int = seq_len
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
                    str(_contig_seq_rec.seq).lower() + \
                    max_num_paddings * PADDING_CHAR
                # note that this does not work for really short contigs
                # e.g. when 'seq_len + paddings > len(contig)'. in this case,
                # the number of paddings will exceed the given maximum
                if len(_padded_seq) < self._seq_len:
                    _warning_msg = \
                        f'The length of contig {_contig_seq_rec.id} from ' \
                        f'genome {_genome_id} is {len(_padded_seq)} with ' \
                        f'paddings at the end, smaller compared to the ' \
                        f'intended dataset sequence length {self._seq_len}' \
                        f'. Ignoring the contig ...'
                    _LOGGER.warning(_warning_msg)
                self._genome_contig_seq_dict[_genome_contig_id] = \
                    _padded_seq
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
        LongTensor and a mask tensor of the same size for padding in boolean
        :rtype: tuple of two tensors
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
                break

        # convert the nucleotide sequence to the indexed (numeric) sequence
        _indexed_seq = torch.LongTensor(
            list(map(NUCLEOTIDE_CHAR_INDEX_DICT.get, _seq)))

        # get the paddings in seq as tensor of booleans
        _padding_mask = (_indexed_seq == PADDING_INDEX)

        return _indexed_seq, _padding_mask
