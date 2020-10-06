"""
File Name:          test_masked_genome_dataset.py
Project:            bioseq-learning

File Description:

"""
import os
import time
import random
import logging
import unittest
from typing import List

from torch.utils.data import DataLoader

from src import E_COLI_GENOME_PARENT_DIR_PATH
from src.datasets import GenomeDataset, GenomeIterDataset
from src.datasets.genome_dataset import \
    PADDING_CHAR, NUCLEOTIDE_CHAR_INDEX_DICT


_PADDING_INDEX: int = NUCLEOTIDE_CHAR_INDEX_DICT[PADDING_CHAR]
_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
_TEST_GENOME_DIR_PATHS: List[str] = [
    os.path.join(
        _TEST_EXAMPLES_DIR_PATH,
        'genome_processing_results_reference',
        '562.2283'
    )]
_TEST_SEQ_LEN: int = 1000
_TEST_NUM_MASKS: int = 5
_TEST_MAX_NUM_PADDINGS: int = 500
_TEST_BATCH_SIZE: int = 16
_TEST_NUM_SAMPLES: int = 10000
_TEST_GENOME_DATASET_KWARGS = {
    'genome_dir_paths': _TEST_GENOME_DIR_PATHS,
    'seq_len': _TEST_SEQ_LEN,
    'max_num_paddings': _TEST_MAX_NUM_PADDINGS,
}

class TestGenomeDataset(unittest.TestCase):
    """unittest class for 'genome_dataset' and 'genome_iter_dataset' classes
    """
    def test_genome_dataset(self):
        """test 'genome_dataset' class
        """
        genome_dataset = GenomeDataset(**_TEST_GENOME_DATASET_KWARGS)
        for _i in range(len(genome_dataset)):
            _indexed_seq, _padding_mask = genome_dataset[_i]
            # test if all the indexed values are in the index dict
            assert set([_v.item() for _v in _indexed_seq.unique()]).issubset(
                set(NUCLEOTIDE_CHAR_INDEX_DICT.values()))
            # test if the number of paddings is valid (no bigger)
            assert _padding_mask.sum().item() <= _TEST_MAX_NUM_PADDINGS

    def test_genome_iter_dataset(self):
        """test 'genome_iter_dataset' class
        """
        # test the strictly sampling dataloader
        # make sure that the number of batches during an iteration only
        # traverse all the samples once, and therefore the total number of
        # batches should be the same as the length of the dataloader
        genome_strict_iter_dataloader = DataLoader(
            GenomeIterDataset(
                **_TEST_GENOME_DATASET_KWARGS,
                strict_iteration=True,
            ),
            batch_size=_TEST_BATCH_SIZE,
            drop_last=True,
        )
        _num_batches: int = 0
        for _batch_data in genome_strict_iter_dataloader:
            assert _batch_data[0].shape == (_TEST_BATCH_SIZE, _TEST_SEQ_LEN)
            assert _batch_data[1].shape == (_TEST_BATCH_SIZE, _TEST_SEQ_LEN)
            _num_batches += 1
            assert _num_batches <= len(genome_strict_iter_dataloader)

        # test the randomly sampling dataloader
        # make sure that the random sampling goes beyond the length of the
        # dataloader (can go on and on but with replacement)
        genome_random_iter_dataloader = DataLoader(
            GenomeIterDataset(**_TEST_GENOME_DATASET_KWARGS),
            batch_size=_TEST_BATCH_SIZE,
            drop_last=True,
        )
        _num_rand_batches: int = 0
        for _batch_data in genome_random_iter_dataloader:
            assert _batch_data[0].shape == (_TEST_BATCH_SIZE, _TEST_SEQ_LEN)
            assert _batch_data[1].shape == (_TEST_BATCH_SIZE, _TEST_SEQ_LEN)
            _num_rand_batches += 1

            if _num_rand_batches > _num_batches:
                return
        assert False

    def test_genome_dataset_indexing_time(self):

        _genome_parent_dir_path: str = E_COLI_GENOME_PARENT_DIR_PATH
        _genome_dir_paths: List[str] = [
            os.path.join(_genome_parent_dir_path, _genome_id)
            for _genome_id in os.listdir(_genome_parent_dir_path)
            if os.path.isdir(os.path.join(_genome_parent_dir_path, _genome_id))
        ]
        logging.getLogger('src.datasets').setLevel(logging.ERROR)
        _genome_dataset = GenomeDataset(
            _genome_dir_paths,
            seq_len=_TEST_SEQ_LEN,
            max_num_paddings=_TEST_MAX_NUM_PADDINGS,
        )

        _test_indices: List[int] = \
            random.sample(range(len(_genome_dataset)), _TEST_NUM_SAMPLES)
        _start_time = time.time()
        for _i in _test_indices:
            _genome_dataset[_i]
        print(f'Indexing {_TEST_NUM_SAMPLES} samples from the dataset '
              f'takes {time.time() - _start_time:.2f} seconds.')


if __name__ == '__main__':
    unittest.main()
