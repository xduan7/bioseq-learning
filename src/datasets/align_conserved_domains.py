"""
File Name:          align_conserved_domains.py
Project:            bioseq-learning

File Description:

"""
import os
import argparse
from multiprocessing import Pool
from typing import List, Tuple, Dict
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from skbio.sequence import Protein
from skbio.alignment import local_pairwise_align_ssw, TabularMSA


BLOSUM62_URL = 'https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt'


def align_conserved_domain_pair(
        cd1_sequence_path: str,
        cd2_sequence_path: str,
        substitution_matrix: dict,
) -> Tuple[str, str, float]:
    """TODO: docstring here

    :param cd1_sequence_path:
    :type cd1_sequence_path:
    :param cd2_sequence_path:
    :type cd2_sequence_path:
    :param substitution_matrix:
    :type substitution_matrix:
    :return:
    :rtype:
    """

    cd1_seq_records = list(SeqIO.parse(cd1_sequence_path, 'fasta'))
    cd2_seq_records = list(SeqIO.parse(cd2_sequence_path, 'fasta'))
    cd1_consensus_seq = str(cd1_seq_records[0].seq)
    cd2_consensus_seq = str(cd2_seq_records[0].seq)

    # Smith Waterman alignment for the consensus sequences
    tabular_msa, alignment_score, position_list = \
        local_pairwise_align_ssw(
        sequence1=Protein(cd1_consensus_seq),
        sequence2=Protein(cd2_consensus_seq),
        substitution_matrix=substitution_matrix,
        protein=True,
    )
    tabular_msa: TabularMSA

    # perform N-to-M alignment and normalize the score with alignment length
    cd1_seq_start, cd1_seq_end = position_list[0][0], position_list[0][1] + 1
    cd2_seq_start, cd2_seq_end = position_list[1][0], position_list[1][1] + 1
    scores = []
    for _cd1_seq_record in cd1_seq_records[1:]:
        for _cd2_seq_record in cd2_seq_records[1:]:
            _, _score, _ = local_pairwise_align_ssw(
                sequence1=Protein(
                    str(_cd1_seq_record.seq)[cd1_seq_start:cd1_seq_end]),
                sequence2=Protein(
                    str(_cd2_seq_record.seq)[cd2_seq_start:cd2_seq_end]),
                substitution_matrix=substitution_matrix,
                protein=True,
            )
            scores.append(_score)

    # return alignment score normalized with alignment length, along with
    # the conserved domain names in a tuple
    cd1_name: str = os.path.basename(cd1_sequence_path).\
        rstrip('.FASTA').rstrip('.fasta')
    cd2_name: str = os.path.basename(cd2_sequence_path).\
        rstrip('.FASTA').rstrip('.fasta')
    return cd1_name, cd2_name, np.average(scores) / tabular_msa.shape[1]


def __align_conserved_domain_pair(arg):
    """wrapper function for align_conserved_domain_pair that only takes
    one tuple argument

    :param arg: tuple argument for align_conserved_domain_pair to be unpacked
    :type arg: Tuple[str, str, dict]
    :return: None
    """

    return align_conserved_domain_pair(*arg)


def align_conserved_domains(
        cd_fasta_parent_dir_path: str,
        alignment_csv_path: str,
        num_workers: int = 1,
):

    # clamp the number of workers between range [1, number of CPU cores]
    num_workers: int = max(1, min(num_workers, os.cpu_count()))

    # prepare for the dictionary substitution matrix (BLOSUM62)
    blosum62_df = pd.read_csv(
        'https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt',
        sep='\s+',
        index_col=0,
        header=0,
        skiprows=6,
    )
    # add a column and a row for gap
    blosum62_df['-'] = 0
    blosum62_df.loc['-'] = 0
    blosum62_dict = blosum62_df.to_dict()

    # get all the paths to FASTA of conserved domains
    cd_fasta_paths: List[str] = []
    for _cd_fasta in os.listdir(cd_fasta_parent_dir_path):
        if _cd_fasta_path.lower().endswith('fasta'):
            _cd_fasta_path: str = \
                os.path.join(cd_fasta_parent_dir_path, _cd_fasta)
            cd_fasta_paths.append(_cd_fasta_path)

    # create a tuple of arguments for align_conserved_domain_pair
    align_conserved_domain_pair_arguments: List[Tuple[str, str, dict]] = \
        [(_paths[0], _paths[1], blosum62_dict)
         for _paths in combinations_with_replacement(cd_fasta_paths, 2)]

    # parallel genome processing with pool
    alignment_results: Dict[str: Dict[str, float]] = {}
    with Pool(num_workers) as _pool:
        for _cd1_name, _cd2_name, _score in tqdm(_pool.imap_unordered(
                __align_conserved_domain_pair,
                align_conserved_domain_pair_arguments),
                total=len(align_conserved_domain_pair_arguments)):
            alignment_results[_cd1_name][_cd2_name] = _score
            alignment_results[_cd2_name][_cd1_name] = _score

    # TODO: write the nested dict to a CSV file with csv.DictWriter
    print(alignment_results)

    # TODO: implement unit test for conserved domain alignment


# script (executable) for conserved domain alignment
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='align conserved domain FASTAs in parallel')

    parser.add_argument(
        '-i', '--cd_fasta_parent_dir_path', type=str,
        help='path to parent directory of the FASTAs of conserved domains)',
    )
    parser.add_argument(
        '-o', '--alignment_csv_path', type=str, default=None,
        help='path to the alignment CSV file',
    )
    parser.add_argument(
        '-t', '--num_workers', type=int, default=1,
        help='number of workers for parallelization',
    )
    args = parser.parse_args()

    # usage example
    # $export PYTHONPATH=<project dir>:$PYTHONPATH
    # $python process_genomes.py \
    #      -i ~/data/PATRIC/genomes \
    #      -o ~/data/PATRIC/processed_genomes \
    #      -t 80
    align_conserved_domains(
        cd_fasta_parent_dir_path=args.cd_fasta_parent_dir_path,
        alignment_csv_path=args.alignment_csv_path,
        num_workers=args.num_workers,
    )
