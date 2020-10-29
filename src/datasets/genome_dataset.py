"""
File Name:          genome_dataset.py
Project:            bioseq-learning-cd

File Description:

"""
import os
import random
import logging
from bisect import bisect
from multiprocessing import Pool, cpu_count
from typing import Tuple, List, Set, Dict, Iterable, Iterator, Optional, Union

import torch
import numpy as np
from Bio import SeqIO, SeqRecord
from torch.utils.data import Dataset, IterableDataset


# MASKED_CHAR: str = '?'
PADDING_CHAR: str = '-'
NUCLEOTIDE_CHAR_INDEX_DICT: Dict[str, int] = {
    PADDING_CHAR: 0,
    'a': 1,
    't': 2,
    'g': 3,
    'c': 4,
    # MASKED_CHAR: 5,
}
INDEX_NUCLEOTIDE_CHAR_DICT: Dict[int, str] = {
    _index: _char for _char, _index in NUCLEOTIDE_CHAR_INDEX_DICT.items()
}
NUCLEOTIDE_CHAR_SET: Set[str] = set(NUCLEOTIDE_CHAR_INDEX_DICT.keys())
NUCLEOTIDE_CHAR_VOCAB_SIZE: int = len(NUCLEOTIDE_CHAR_INDEX_DICT)
PADDING_INDEX: int = NUCLEOTIDE_CHAR_INDEX_DICT[PADDING_CHAR]

_LOGGER = logging.getLogger(__name__)


def _process_single_contig(
        _genome_id: str,
        _contig_seq_path: str,
        _seq_len: int,
        _max_num_paddings: int,
        _nucleotide_char_index_dict: Dict[str, int],
        _padding_index: int,
) -> Optional[Tuple[str, int, np.ndarray, np.ndarray]]:
    """process a single contig sequence, which includes sanity check for
    length and characters, and then transform/index the contig into indices
    that can be digest by word embedding layers

    :param _genome_id: ID of the genome (NOT contig); the contig ID will be
    parsed and appended to genome ID for the ultimate genome contig ID
    :type _genome_id: str
    :param _contig_seq_path: path to the contig sequence fasta file
    :type _contig_seq_path: str
    :param _seq_len: desired length of the output segmented sequence; used
    for sanity checking
    :type _seq_len: int
    :param _max_num_paddings: maximum number of paddings for a output
    sequence; used for sanity checking and sequence construction
    :type _max_num_paddings: int
    :param _nucleotide_char_index_dict: dictionary that maps nucleotide
    characters (ATGCs) to indices
    :type _nucleotide_char_index_dict: Dict[str, np.uint8]
    :param _padding_index: index for the padding character
    :type _padding_index: np.uint8
    :return: if contig sequence is valid in terms of length and nucleotide
    bases, then return a tuple of ID, number of sequences for the contig in
    total, and numpy array of indexes sequence and the padding mask;
    otherwise return None
    :rtype: Optional[Tuple[str, int, np.ndarray, np.ndarray]]
    """

    with open(_contig_seq_path, 'r') as _fh:
        _contig_seq_rec: SeqRecord = \
            next(SeqIO.parse(_fh, 'fasta'))
    _genome_contig_id: str = \
        f'{_genome_id}/{_contig_seq_rec.id}'
    _seq: str = str(_contig_seq_rec.seq).lower()

    # note that this does not work for really short contigs
    # e.g. when 'seq_len + paddings > len(contig)'. in this case,
    # the number of paddings will exceed the given maximum
    if len(_seq) < _seq_len - _max_num_paddings:
        _warning_msg = \
            f'The length of contig {_contig_seq_rec.id} from genome ' \
            f'{_genome_id} is {len(_seq)} with paddings at the end, ' \
            f'smaller compared to the intended dataset sequence ' \
            f'length {_seq_len}. Ignoring the contig ...'
        _LOGGER.warning(_warning_msg)
        return None
    else:
        _num_seqs: int = len(_seq) + _max_num_paddings - _seq_len + 1

    # note that there are sequences with letters are not atgc
    # e.g. genome/contig ID: 562.54941/CP047077
    try:
        _indexed_seq: np.ndarray = np.array(
            [_nucleotide_char_index_dict[_c] for _c in _seq] +
            [_padding_index] * _max_num_paddings,
            dtype=np.uint8
        )
        _padding_mask_seq: np.ndarray = np.array(
            [False] * len(_seq) + [True] * _max_num_paddings,
            dtype=np.bool
        )
        # transforming the indexed sequence into PyTorch tensor might
        # cause: 'RuntimeError: received 0 items of ancdata'
        # _indexed_seq_tensor: torch.LongTensor = \
        #     torch.LongTensor(_indexed_seq)
    except KeyError:
        _warning_msg = \
            f'The contig {_contig_seq_rec.id} from genome {_genome_id} ' \
            f'contains the following letters that are not normal ' \
            f'nucleotide bases: {set(_seq) - NUCLEOTIDE_CHAR_SET}' \
            f'. Ignoring the contig ...'
        _LOGGER.warning(_warning_msg)
        return None
    return _genome_contig_id, _num_seqs, _indexed_seq, _padding_mask_seq


class GenomeDataset(Dataset):
    """basic genome dataset for language model
    """
    def __init__(
            self,
            genome_dir_paths: Iterable[str],
            seq_len: int,
            max_num_paddings: int,
            dry_run: bool,
            dry_run_contig_len: int,
            dry_run_num_contigs: int,
            dry_run_num_ins_per_contig: int,
    ):
        """create a simple dataset of segmented genome sequences

        The dataset will parse all the contigs under all given paths to
        genomes. After checking if the contigs are legal (only includes ATGCs
        for now), the dataset will cut them into given length with overlaps
        and paddings. For example, a sequence of 'atgcatgc' will be
        segmented into 'atgca', tgcat', 'gcatg', 'catgc','atgc-', 'tgc--'
        when seq_len=5 and max_num_paddings=2.


        :param genome_dir_paths: an iterable of paths to individual genomes
        that are processed (check process_genome.py for details)
        :type genome_dir_paths: iterable of strings
        :param seq_len: the length of genome sequences of the dataset after
        segmentation and padding
        :type seq_len: int
        :param max_num_paddings: maximum number of padding characters in a
        genome sequence
        :type max_num_paddings: int
        :param dry_run: TODO
        :type dry_run: bool
        :param dry_run_contig_len: TODO
        :type dry_run_contig_len: int
        :param dry_run_num_contigs: TODO
        :type dry_run_num_contigs: int
        :param dry_run_num_ins_per_contig: TODO
        :type dry_run_num_ins_per_contig: int
        """
        # sanity check for the sequence length and number of paddings
        assert 0 <= max_num_paddings < seq_len

        self._len: int = 0
        self._seq_len: int = seq_len
        self._max_num_paddings: int = max_num_paddings

        # get a list of argument for single contig processing
        # _process_single_contig_arg_list: \
        #     List[Tuple[str, str, int, int, Dict[str, int], int]] = []
        _process_single_contig_arg_list: \
            List[Tuple[str, str, int, int, Dict[str, int], int]] = []
        for _genome_dir_path in genome_dir_paths:

            _genome_id: str = os.path.basename(_genome_dir_path.rstrip('/'))
            _genome_contig_seq_dir_path: str = os.path.join(
                _genome_dir_path, 'contigs')

            for _contig_seq_file_name in \
                    os.listdir(_genome_contig_seq_dir_path):

                _contig_seq_path = os.path.join(
                    _genome_contig_seq_dir_path,
                    _contig_seq_file_name,
                )
                if os.path.isfile(_contig_seq_path):
                    _process_single_contig_arg_list.append((
                        _genome_id,
                        _contig_seq_path,
                        self._seq_len,
                        self._max_num_paddings,
                        # a copy of dict can boost the performance,
                        # probably due to faster memory access for
                        # all the processes during multiprocessing
                        NUCLEOTIDE_CHAR_INDEX_DICT.copy(),
                        PADDING_INDEX,
                    ))

        # multi-processing the genome contigs
        with Pool(cpu_count()) as _pool:
            _processed_contigs: \
                List[Optional[Tuple[str, int, np.ndarray, np.ndarray]]] = \
                _pool.starmap(
                    _process_single_contig,
                    _process_single_contig_arg_list,
                )

        # dict that maps (genome contig id) to the num of sequences that can
        # be segmented from the whole contig, the indexes sequence, and the
        # padding mask for the whole contig sequence
        self._genome_contig_seq_dict: \
            Dict[str, Tuple[int, np.ndarray, np.ndarray]] = {}

        # TODO: debug or info for how many contigs are kept and thrown

        # pick a single contig with about the certain length
        if dry_run:
            # TODO: put this part in another function
            _dry_run_contig: Tuple[str, int, np.ndarray, np.ndarray]
            # _min_contig_len_diff: Union[int, float] = float('inf')
            _max_contig_len: int = 0
            for _processed_single_contig in _processed_contigs:
                if not _processed_single_contig:
                    continue

                _, _, _indexed_seq, _padding_mask_seq \
                    = _processed_single_contig
                _genome_contig_len = \
                    len(_indexed_seq) - self._max_num_paddings
                # if abs(dry_run_contig_len - _genome_contig_len) \
                #         < _min_contig_len_diff:
                if _genome_contig_len > _max_contig_len:
                    # _min_contig_len_diff = \
                    #     abs(dry_run_contig_len - _genome_contig_len)
                    _max_contig_len = _genome_contig_len
                    _dry_run_contig = _processed_single_contig

            for __i in range(dry_run_num_contigs):

                __s = dry_run_contig_len * __i + __i
                __e = dry_run_contig_len * (__i + 1) + __i \
                    - dry_run_num_ins_per_contig

                __k = f'{_dry_run_contig[0]}({__i + 1:2d})'
                __indexed_seq = _dry_run_contig[2][__s:__e]
                __padding_mask_seq = _dry_run_contig[3][__s:__e]

                for _ in range(dry_run_num_ins_per_contig):
                    __k = random.randint(0, len(__indexed_seq))
                    __n = random.sample([1, 2, 3, 4], 1)[0]
                    __indexed_seq = np.insert(__indexed_seq, __k, __n)
                    __padding_mask_seq \
                        = np.insert(__padding_mask_seq, __k, False)

                self._genome_contig_seq_dict[__k] = (
                    dry_run_contig_len - self._seq_len + 1,
                    __indexed_seq,
                    __padding_mask_seq,
                )

        else:
            for _processed_single_contig in _processed_contigs:
                # returned None from process single contig means either length
                # too short or contains illegal character
                if not _processed_single_contig:
                    continue

                # otherwise add the genome contigs into class dict
                _genome_contig_id, _num_seq, _indexed_seq, _padding_mask_seq \
                    = _processed_single_contig
                self._genome_contig_seq_dict[_genome_contig_id] = \
                    (_num_seq, _indexed_seq, _padding_mask_seq)

        # get the tuple of genome contig IDs, and the accumulative number of
        # sequences of all the genome contigs; for example:
        # genome_contig_tuple = (
        #     'contig_0',
        #     'contig_1', ...
        # )
        # acc_num_seqs_tuple = (
        #     num seq of 'contig_0',
        #     num seq of 'contig_0' + num seq of 'contig_1', ...
        # )
        _total_num_seqs: int = 0
        _acc_num_seqs_list: List[int] = []
        _genome_contig_id_list: List[str] = []
        for _genome_contig_id, _seq_tuple in \
                self._genome_contig_seq_dict.items():
            _num_seqs = _seq_tuple[0]
            _total_num_seqs += _num_seqs
            _acc_num_seqs_list.append(_total_num_seqs)
            _genome_contig_id_list.append(_genome_contig_id)
        self._acc_num_seqs_tuple: Tuple[int, ...] = \
            tuple(_acc_num_seqs_list)
        self._genome_contig_id_tuple: Tuple[str, ...] = \
            tuple(_genome_contig_id_list)

        self._len = _total_num_seqs

    def __len__(self) -> int:
        """get the length of the dataset

        :return: number of samples in total
        :rtype: int
        """
        return self._len

    def __getitem__(
            self,
            index: int,
    ) -> Tuple[torch.LongTensor, torch.BoolTensor]:
        """
        :param index: index to the genome sequence in [0, len(self)]
        :type index: int
        :return: a tuple that contains a indexed genome sequence in
        LongTensor and a mask tensor of the same size for padding in boolean
        :rtype: tuple of two tensors
        """
        # for negative indexing that is probably not going to be used ...
        _index: int = (index + self._len) if index < 0 else index

        # bisect searching reduce the time from 7.92 seconds per 10000
        # samples to 0.45 seconds
        # get the index for both tuples genome contig IDs and the
        # accumulative numbers of the sequences, with bisect
        _genome_contig_index: int = \
            bisect(self._acc_num_seqs_tuple, _index)
        # truncate the index to get the starting position of the sequence,
        # which is then found with the tuple of genome contig IDs
        # then the whole genome sequence, and finally segmented into the
        # actual indexed sequence that we are looking for ...
        _seq_start_pos: int = _index if _genome_contig_index == 0 else \
            (_index - self._acc_num_seqs_tuple[_genome_contig_index - 1])
        _seq_end_pos: int = _seq_start_pos + self._seq_len
        _genome_contig_id: str = \
            self._genome_contig_id_tuple[_genome_contig_index]

        # segment the indexed sequence and padding mask sequence
        _indexed_seq: np.ndarray
        _padding_mask_seq: np.ndarray
        _, _indexed_seq, _padding_mask_seq = \
            self._genome_contig_seq_dict[_genome_contig_id]

        _segmented_indexed_seq: np.ndarray = \
            _indexed_seq[_seq_start_pos: _seq_end_pos]
        _segmented_indexed_seq_tensor: torch.LongTensor = \
            torch.LongTensor(_segmented_indexed_seq)

        _segmented_padding_mask_seq: np.ndarray = \
            _padding_mask_seq[_seq_start_pos: _seq_end_pos]
        _segmented_padding_mask_tensor: torch.BoolTensor = \
            torch.BoolTensor(_segmented_padding_mask_seq)

        return _segmented_indexed_seq_tensor, _segmented_padding_mask_tensor


class GenomeIterDataset(IterableDataset, GenomeDataset):
    """basic iterative genome dataset for language model
    """

    def __init__(
            self,
            genome_dir_paths: Iterable[str],
            seq_len: int,
            max_num_paddings: int,
            strict_iteration: bool,
            dry_run: bool,
            dry_run_contig_len: int,
            dry_run_num_contigs: int,
            dry_run_num_ins_per_contig: int,
    ):
        """create a genome iterable dataset of segmented genome sequences

        This dataset is pretty much the same as the base GenomeDataset,
        except that it's iterable. When a PyTorch dataloader takes such
        iterable dataset, it won't reach all the way to the end of the
        dataset. Instead, each batch is generated with the __next__ method
        on the fly. This could save us a ton of time during the
        initialization of dataloader.

        :param genome_dir_paths: an iterable of paths to individual genomes
        that are processed (check process_genome.py for details)
        :type genome_dir_paths: iterable of strings
        :param seq_len: the length of genome sequences of the dataset after
        segmentation and padding
        :type seq_len: int
        :param max_num_paddings: maximum number of padding characters in a
        genome sequence
        :type max_num_paddings: int
        :param strict_iteration: indicator for strict iteration, which means
        that for a single iteration, an index will not be traversed twice;
        if set to False, the iterative dataset will randomly grab a sample
        without any bookkeeping action
        :type strict_iteration: bool
        :param dry_run: TODO
        :type dry_run: bool
        :param dry_run_contig_len: TODO
        :type dry_run_contig_len: int
        :param dry_run_num_contigs: TODO
        :type dry_run_num_contigs: int
        :param dry_run_num_ins_per_contig: TODO
        :type dry_run_num_ins_per_contig: int
        """
        super(GenomeIterDataset, self).__init__(
            genome_dir_paths=genome_dir_paths,
            seq_len=seq_len,
            max_num_paddings=max_num_paddings,
            dry_run=dry_run,
            dry_run_contig_len=dry_run_contig_len,
            dry_run_num_contigs=dry_run_num_contigs,
            dry_run_num_ins_per_contig=dry_run_num_ins_per_contig,
        )
        self._strict_iteration: bool = strict_iteration
        if self._strict_iteration:
            _warning_msg = \
                f'PyTorch IterableDataset with strict iteration is ' \
                f'functionally the same as Dataset, but slower due to the ' \
                f'poor handling of shuffling. Use GenomeDataset and ' \
                f'shuffle=True for its DataLoader if necessary.'
            _LOGGER.warning(_warning_msg)

    def __refresh_indices(self) -> None:
        """refresh the shuffled indices for a new iteration
        """
        self._curr_pos: int = 0
        self._shuffled_indices: List[int] = list(range(self._len))
        random.shuffle(self._shuffled_indices)

    def __strictly_next(self) -> Tuple[torch.LongTensor, torch.BoolTensor]:
        """get the next sample from the dataset, which is indexed with a
        shuffled list and a position pointer

        :return: a tuple that contains a indexed genome sequence in
        LongTensor and a mask tensor of the same size for padding in boolean
        :rtype: tuple of two tensors
        """
        try:
            _curr_index: int = self._shuffled_indices[self._curr_pos]
            self._curr_pos = self._curr_pos + 1
            return self[_curr_index]
        except IndexError:
            _error_msg = f'reached the end of the genome dataset'
            raise StopIteration(_error_msg)

    def __randomly_next(self) -> Tuple[torch.LongTensor, torch.BoolTensor]:
        """get the a random sample from the dataset

        :return: a tuple that contains a indexed genome sequence in
        LongTensor and a mask tensor of the same size for padding in boolean
        :rtype: tuple of two tensors
        """
        return self[random.randint(0, len(self) - 1)]

    def __iter__(self) -> Iterator:
        """initialize an iterator for the dataset

        :return: iterable of the dataset
        :rtype: Iterator
        """
        if self._strict_iteration:
            self.__refresh_indices()
        return self

    def __next__(self) -> Tuple[torch.LongTensor, torch.BoolTensor]:
        """get the next sample from the dataset

        :return: a tuple that contains a indexed genome sequence in
        LongTensor and a mask tensor of the same size for padding in boolean
        :rtype: tuple of two tensors
        """
        if self._strict_iteration:
            return self.__strictly_next()
        else:
            return self.__randomly_next()
