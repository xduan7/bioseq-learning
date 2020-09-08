"""
File Name:          test_process_genomes.py
Project:            bioseq-learning

File Description:

"""
import os
import unittest
from io import StringIO
from filecmp import dircmp

from src.datasets import process_genomes


_TEST_EXAMPLES_DIR_PATH: str = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
_GENOME_PARENT_DIR_PATH: str = \
    os.path.join(_TEST_EXAMPLES_DIR_PATH, 'genomes')
_OUTPUT_PARENT_DIR_PATH: str = \
    os.path.join(_TEST_EXAMPLES_DIR_PATH, 'genome_processing_results')
_REFERENCE_OUTPUT_PARENT_DIR_PATH: str = \
    os.path.join(
        _TEST_EXAMPLES_DIR_PATH,
        'genome_processing_results_reference',
    )


class TestProcessGenomes(unittest.TestCase):
    """unittest class for 'process_genomes' function
    """
    def test_process_genomes(self):
        """test 'process_genomes' function
        """
        process_genomes(_GENOME_PARENT_DIR_PATH, _OUTPUT_PARENT_DIR_PATH, 4)

        # compare the output with the reference
        from contextlib import redirect_stdout
        _str_io: StringIO = StringIO()
        with redirect_stdout(_str_io):
            dircmp(
                _OUTPUT_PARENT_DIR_PATH,
                _REFERENCE_OUTPUT_PARENT_DIR_PATH,
            ).report_full_closure()
        assert 'Differing files' not in _str_io.getvalue()


if __name__ == '__main__':
    unittest.main()
