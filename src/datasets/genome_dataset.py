"""
File Name:          genome_dataset.py
Project:            bioseq-learning-cd

File Description:

"""
from typing import Union
from torch.utils.data import Dataset


class BasicGenomeMLMDataset(Dataset):
    """
    Basic genome dataset for (dynamic) masked language model training

    """

    def __init__(
            self,
            genome_dir: str,
            seq_len: int,
            num_masks: Union[int, float],
            max_contigs: int = 32,
            # minimum number of real nucleotide?
            # indicator and function of weighing masks?
    ):
        # TODO: search for standard NLP example for MLM
        # (1) walk through genome directory
        # (2) get the total number of test_process_genomes
        # (3) get the number of contigs and lengths for each genome
        # (4) calculate the total number of data sequences
        # (5) index the data sequences, but how?
        pass

    def __len__(self):
        pass

    def __getitem__(
            self,
            index: int,
    ):
        # apply mask dynamically
        # get the weights of the masks if required, based on the genomic
        # feature annotations and conserved domain search
        pass
