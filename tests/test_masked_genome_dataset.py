"""
File Name:          test_masked_genome_dataset.py
Project:            bioseq-learning

File Description:

"""
import os
import unittest
from typing import List, Union

import numpy as np

from src.datasets import MaskedGenomeDataset
from src.datasets.masked_genome_dataset import \
    PADDING_CHAR, NUCLEOTIDE_CHAR_INDEX_DICT


_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
_TEST_GENOME_DIR_PATHS: List[str] = [
    os.path.join(
        _TEST_EXAMPLES_DIR_PATH,
        'genome_processing_results_reference',
        '562.2283'
    )]
_TEST_SEQ_LEN: int = 100
_TEST_NUM_MASKS: int = 5
_TEST_MAX_NUM_PADDINGS: int = 10
_PADDING_INDEX: int = NUCLEOTIDE_CHAR_INDEX_DICT[PADDING_CHAR]


class TestMaskedGenomeDataset(unittest.TestCase):
    """unittest class for 'masked_genome_dataset' class
    """
    def test_masked_genome_dataset(self):
        """test 'masked_genome_dataset' class
        """
        masked_genome_dataset = MaskedGenomeDataset(
            genome_dir_paths=_TEST_GENOME_DIR_PATHS,
            seq_len=_TEST_SEQ_LEN,
            num_masks=_TEST_NUM_MASKS,
            max_num_paddings=_TEST_MAX_NUM_PADDINGS,
        )

        _indexed_seq, _mask = masked_genome_dataset[0]

        # make sure that the number of masks (non-masks) is correct
        _mask_unique_values, _mask_unique_value_counts = \
            _mask.unique(return_counts=True)
        _non_mask_unique_value_index = \
            np.argwhere(_mask_unique_values.numpy() == 0.).flatten()[0]
        assert _mask_unique_value_counts[_non_mask_unique_value_index].item() \
               == (_TEST_SEQ_LEN - _TEST_NUM_MASKS)

        # make sure that the paddings are correct
        for _i in range(_TEST_MAX_NUM_PADDINGS):
            assert _indexed_seq[_i].item() == _PADDING_INDEX

        _indexed_seq, _mask = masked_genome_dataset[-1]

        # make sure that the paddings are correct when paddings are at the end
        for _i in range(_TEST_MAX_NUM_PADDINGS):
            assert _indexed_seq[-(_i + 1)].item() == _PADDING_INDEX


if __name__ == '__main__':
    unittest.main()
