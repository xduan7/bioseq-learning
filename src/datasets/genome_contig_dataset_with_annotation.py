"""
File Name:          genome_dataset_with_annotation.py
Project:            bioseq-learning

File Description:

"""
import os
import random
import logging
import resource

from typing import Dict, Iterable, Optional, Sequence, Tuple

import torch
from torch.utils.data import Dataset, IterableDataset


_LOGGER = logging.getLogger(__name__)
# prevents the OS from raising "too many files open" error during the
# parallelized loading and processing of the sequences
resource.setrlimit(
    resource.RLIMIT_NOFILE,
    (4096, resource.getrlimit(resource.RLIMIT_NOFILE)[1]),
)


class GenomeContigDataset(Dataset):

    def __init__(
            self,
            genome_contig_paths: Iterable[str],
            seq_len: int,
            kmer_len: int,
            include_cr_info: bool,
            include_cd_info: bool,
            max_num_contigs: Optional[int] = None,
            max_num_paddings: Optional[int] = None,
            contig_len_range:
            Tuple[Optional[int], Optional[int]] = (None,  None),
    ):
        """There are two steps from the paths of genome contigs to the 
        inputs/outputs tensors of neural networks:
        (1) (during init) genome contig path - > genome contig seq;
        (2) (during indexing) genome contig seq -> indexed kmer masked seq;
        
        Note that the indexing might lead to a huge load for CPU, 
        which might cause lower utilization for GPUs.
        
        :param genome_contig_paths: 
        :type genome_contig_paths: 
        :param seq_len: 
        :type seq_len: 
        :param kmer_len: 
        :type kmer_len: 
        :param include_cr_info: 
        :type include_cr_info: 
        :param include_cd_info: 
        :type include_cd_info: 
        :param max_num_contigs: 
        :type max_num_contigs: 
        :param max_num_paddings: 
        :type max_num_paddings: 
        :param contig_len_range: 
        :type contig_len_range: 
        """

        # number of indices with given kmer (redundant?)
        self.num_indices: int

        # dictionaries that map genome sequence to k-mer sequence and back
        # TODO: callables that uses the dictionaries and maps seq <-> kmer seq
        self.seq_to_kmer_seq_dict: Dict[Sequence[int], int]
        self.kmer_seq_to_seq_dict: Dict[int, Sequence[int]]

        pass

    def __getitem__(
            self,
            index: int,
    ) -> Tuple[torch.LongTensor, torch.LongTensor, torch.BoolTensor]:
        # return the tuple of
        # indexed_kmer_masked_seq <-- input
        # indexed_kmer_seq <-- target
        # boolean seq for coding region (CR) annotation
        # TODO: conserved domain (CD) annotation
        pass
