"""
File Name:          test_masked_genome_dataset.py
Project:            bioseq-learning

File Description:

"""
import os
import unittest

from src.datasets import MaskedGenomeDataset


_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
_REFERENCE_OUTPUT_PARENT_DIR_PATH: str = \
    os.path.join(_TEST_EXAMPLES_DIR_PATH,
                 'genome_processing_results_reference')


class TestMaskedGenomeDataset(unittest.TestCase):
    """unittest class for 'masked_genome_dataset' class
    """
    def test_masked_genome_dataset(self):
        """test 'masked_genome_dataset' class
        """
        masked_genome_dataset = MaskedGenomeDataset(
            genome_parent_dir_path=_REFERENCE_OUTPUT_PARENT_DIR_PATH,
            seq_len=100,
            num_masks=0.05,
            max_num_paddings=10,
        )

        print(masked_genome_dataset[0])
        print(masked_genome_dataset[0])
        print(masked_genome_dataset[1])


if __name__ == '__main__':
    unittest.main()
