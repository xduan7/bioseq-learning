"""
File Name:          test_process_genomes.py
Project:            bioseq-learning

File Description:

"""
import unittest


_GENOME_PARENT_DIR_PATH: str = 'examples/genomes/'
_REFERENCE_RESULTS_DIR_PATH: str = \
    'examples/genome_processing_reference_results/'


class TestProcessGenomes(unittest.TestCase):
    """unittest class for 'process_genomes' function
    """
    def test_process_genomes(self):
        """test 'process_genomes' function
        """
        # TODO: add small genome data for testing
        # TODO: test genome processing code here (checking results how)
        pass


if __name__ == '__main__':
    unittest.main()
