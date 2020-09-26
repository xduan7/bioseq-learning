"""
File Name:          test_masked_genome_dataset.py
Project:            bioseq-learning

File Description:

"""
import os
import unittest
from typing import List

from src.datasets import MaskedGenomeDataset


_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
# _REFERENCE_OUTPUT_PARENT_DIR_PATH: str = \
#     os.path.join(_TEST_EXAMPLES_DIR_PATH,
#                  'genome_processing_results_reference')

_TEST_GENOME_DIR_PATHS: List[str] = [
    os.path.join(
        _TEST_EXAMPLES_DIR_PATH,
        'genome_processing_results_reference',
        '562.2283'
    )
]


class TestMaskedGenomeDataset(unittest.TestCase):
    """unittest class for 'masked_genome_dataset' class
    """
    def test_masked_genome_dataset(self):
        """test 'masked_genome_dataset' class
        """
        masked_genome_dataset = MaskedGenomeDataset(
            genome_dir_paths=_TEST_GENOME_DIR_PATHS,
            seq_len=100,
            num_masks=0.05,
            max_num_paddings=10,
        )

        # TODO: write better test cases if necessary
        _indexed_seq, _mask = masked_genome_dataset[0]
        print(_indexed_seq)
        print(_mask)
        _indexed_seq, _mask = masked_genome_dataset[0]
        print(_indexed_seq)
        print(_mask)
        _indexed_seq, _mask = masked_genome_dataset[1]
        print(_indexed_seq)
        print(_mask)


if __name__ == '__main__':
    unittest.main()
