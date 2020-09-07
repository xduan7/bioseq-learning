"""
File Name:          test_process_genomes.py
Project:            bioseq-learning

File Description:

"""
import os
import unittest

from src.datasets import process_genomes


_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
_GENOME_PARENT_DIR_PATH: str = \
    os.path.join(_TEST_EXAMPLES_DIR_PATH, 'genomes')
# _REFERENCE_RESULTS_DIR_PATH: str = \
#     os.path.join(
#         _TEST_EXAMPLES_DIR_PATH,
#         'genome_processing_reference_results',
#     )


class TestProcessGenomes(unittest.TestCase):
    """unittest class for 'process_genomes' function
    """
    def test_process_genomes(self):
        """test 'process_genomes' function
        """
        process_genomes(_GENOME_PARENT_DIR_PATH, 4)


if __name__ == '__main__':
    unittest.main()
