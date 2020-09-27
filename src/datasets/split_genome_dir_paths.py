"""
File Name:          split_genome_dir_paths.py
Project:            bioseq-learning

File Description:

    This file implements a function for genome directory paths
    (training/validation/testing) splitting.

"""
import os
from random import shuffle
from typing import Tuple, List


def split_genome_dir_paths(
        genome_parent_dir_path: str,
        vld_ratio: float,
        tst_ratio: float,
) -> Tuple[List[str], List[str], List[str]]:
    """split the individual genomes paths into three folds, for training,
    validation, and testing separately

    :param genome_parent_dir_path: path to the parent directory of all the
    genomes for the learning task, note that the folder names of each
    individual genome will be used as genome ID, and each folder must
    contain 'contigs', 'features', and 'conserved domains' directories
    :type genome_parent_dir_path: str
    :param vld_ratio: ratio of validation set (out of all genomes)
    :type vld_ratio: float
    :param tst_ratio: ratio of test set (out of all genomes)
    :type tst_ratio: float
    :return: training, validation, and test sets of paths to genomes in order
    :rtype: tuple of three lists of strings
    """

    # get the paths for all the genomes
    genome_dir_paths: List[str] = [
        os.path.join(genome_parent_dir_path, _genome_id)
        for _genome_id in os.listdir(genome_parent_dir_path)
        if os.path.isdir(os.path.join(genome_parent_dir_path, _genome_id))
    ]
    shuffle(genome_dir_paths)

    # get the number of samples for validation and testing
    _num_samples: int = len(genome_dir_paths)
    num_vld_samples: int = round(_num_samples * vld_ratio)
    num_tst_samples: int = round(_num_samples * tst_ratio)
    num_trn_samples: int = _num_samples - num_vld_samples - num_tst_samples
    assert num_vld_samples >= 0
    assert num_tst_samples >= 0
    assert num_trn_samples >= 0

    return (
        genome_dir_paths[:num_trn_samples],
        genome_dir_paths[num_trn_samples:num_trn_samples + num_vld_samples],
        genome_dir_paths[num_trn_samples + num_vld_samples:],
    )
