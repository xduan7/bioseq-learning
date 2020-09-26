"""
File Name:          process_genomes.py
Project:            bioseq-learning-cd

File Description:

"""
import os
import json
import logging
import argparse
from multiprocessing import Pool
from typing import Optional, Iterator, List, Tuple

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO, SeqRecord

from src.utilities import create_directory
from src.datasets import conserved_domain_search


_LOGGER = logging.getLogger(__name__)


def process_genome(
        genome_dir_path: str,
        output_dir_path: Optional[str] = None,
        genome_id: Optional[str] = None,
        no_cd_search: bool = False,
):
    """process a PATRIC genome

    the processing steps are listed as the following:
    - prepare the directories
    - split all the contigs and their features
    - search for conserved domains on all the contigs
    - create an info file

    :param genome_dir_path: path to PATRIC genome directory
    :type genome_dir_path: str
    :param output_dir_path: optional parent path to the processed genomes;
    use the genome directory if not given
    :type output_dir_path: str
    :param genome_id: optional ID for genome ID; using the genome
    directory name if not given
    :type genome_id: str
    :param no_cd_search: indicator for disabling conserved domain search,
    which means parsing the contig and genome information only
    :type no_cd_search: bool
    :return: None
    """

    # get the genome ID if not given
    genome_dir_path: str = os.path.join(genome_dir_path, '')
    genome_id: str = genome_id if genome_id else \
        os.path.basename(genome_dir_path[:-1])

    # get the genome and feature paths
    contig_seq_path: str = os.path.join(genome_dir_path, f'{genome_id}.fna')
    feature_path: str = \
        os.path.join(genome_dir_path,  f'{genome_id}.PATRIC.features.tab')
    if not (os.path.exists(contig_seq_path) and os.path.exists(feature_path)):
        _warning_msg = \
            f'Genome {genome_id} has incomplete data: ' \
            f'missing nucleotide sequence file ({genome_id}.fna) and/or ' \
            f'PATRIC feature table file ({genome_id}.PATRIC.features.tab).'
        _LOGGER.warning(_warning_msg)
        return

    # create the subdirectories if not exist
    # - ./OUTPUT_DIR: output directory if given
    # - ./OUTPUT_DIR/contigs: contig nucleotide sequences
    # - ./OUTPUT_DIR/features:  feature annotations for every contig
    # - ./OUTPUT_DIR/conserved_domains: conserved domains for every contig
    output_dir_path: str = output_dir_path \
        if output_dir_path else genome_dir_path
    contig_dir_path: str = os.path.join(output_dir_path, 'contigs')
    feature_dir_path: str = os.path.join(output_dir_path, 'features')
    conserved_domain_dir_path: str = \
        os.path.join(output_dir_path, 'conserved_domains')
    create_directory(output_dir_path)
    create_directory(contig_dir_path)
    create_directory(feature_dir_path)
    create_directory(conserved_domain_dir_path)

    # read the entire genome into list of contig sequences
    contig_seq_record_iter: Iterator[SeqRecord] = \
        SeqIO.parse(contig_seq_path, 'fasta')

    # split genome contigs and write into files named with their PATRIC IDs
    contig_info_dict: dict = {}
    _contig_seq_record: SeqRecord
    for _contig_seq_record in contig_seq_record_iter:
        _contig_id: str = _contig_seq_record.id
        _contig_seq_path: str = \
            os.path.join(contig_dir_path, f'{_contig_id}.fna')
        contig_info_dict[_contig_id] = {
            'length': len(_contig_seq_record),
        }
        with open(_contig_seq_path, 'w+') as _f:
            SeqIO.write(_contig_seq_record, _f, 'fasta')

    # split PATRIC features and write into files named with their PATRIC IDs
    # assume that accession is the same as contig in PATRIC
    feature_df = pd.read_table(feature_path)
    _accession: str
    _feature_by_accession_df: pd.DataFrame
    for _accession, _feature_by_accession_df \
            in feature_df.groupby('accession'):
        _feature_by_accession_path: str = \
            os.path.join(feature_dir_path, f'{_accession}.tsv')
        _source_mask = (_feature_by_accession_df['feature_type'] == 'source')
        _feature_by_accession_df[~_source_mask].to_csv(
            _feature_by_accession_path, index=None, sep='\t')

    # iterate through all the contigs again for conserved domain search
    if not no_cd_search:
        for _contig_id in contig_info_dict.keys():
            _contig_seq_path: str = \
                os.path.join(contig_dir_path, f'{_contig_id}.fna')
            _contig_cd_ans_path: str = \
                os.path.join(conserved_domain_dir_path, f'{_contig_id}.ans')
            _contig_cd_xml_path: str = \
                os.path.join(conserved_domain_dir_path, f'{_contig_id}.xml')
            _contig_cd_txt_path: str = \
                os.path.join(conserved_domain_dir_path, f'{_contig_id}.txt')
            _contig_cd_csv_path: str = \
                os.path.join(conserved_domain_dir_path, f'{_contig_id}.csv')
            conserved_domain_search(
                nucleotide_seq_path=_contig_seq_path,
                cd_ans_path=_contig_cd_ans_path,
                cd_xml_path=_contig_cd_xml_path,
                cd_txt_path=_contig_cd_txt_path,
                cd_csv_path=_contig_cd_csv_path,
            )

    # save the genome information in json format
    genome_info_path: str = os.path.join(output_dir_path, 'info.json')
    genome_info: dict = {
        'genome_name': feature_df['genome_name'].unique()[0],
        'number_of_contigs': len(contig_info_dict),
        'contigs': contig_info_dict,
    }
    with open(genome_info_path, 'w+') as _fh:
        json.dump(genome_info, _fh, indent=4)


def __process_genome(arg):
    """wrapper function for process_genome that only takes one tuple argument

    :param arg: tuple argument for process_genome to be unpacked
    :type arg: Tuple[str, Optional[str], Optional[str]]
    :return: None
    """

    return process_genome(*arg)


def process_genomes(
        genome_parent_dir_path: str,
        output_parent_dir_path: Optional[str] = None,
        no_cd_search: bool = False,
        num_workers: int = 1,
):
    """process PATRIC genomes in the given parent directory in parallel

    :param genome_parent_dir_path: path to the parent directory of all
    the PATRIC genome directories to be processed
    :type genome_parent_dir_path: str
    :param output_parent_dir_path: optional path to the parent directory of
    all the processed genomes
    :type output_parent_dir_path: str
    :param num_workers: maximum number of workers for parallelization
    :type num_workers: int
    :return: None
    """

    # clamp the number of workers between range [1, number of CPU cores]
    num_workers: int = max(1, min(num_workers, os.cpu_count()))

    # get all the paths to genomes, and output paths if possible
    process_genome_arguments: List[Tuple[str, Optional[str], str, bool]] = []
    for _genome_id in os.listdir(genome_parent_dir_path):
        _genome_dir_path: str = \
            os.path.join(genome_parent_dir_path, _genome_id)
        if os.path.isdir(_genome_dir_path):
            _output_dir_path: Optional[str] = \
                os.path.join(output_parent_dir_path, _genome_id) \
                if output_parent_dir_path else None
            process_genome_arguments.append(
                (_genome_dir_path, _output_dir_path, _genome_id, no_cd_search))

    # parallel genome processing with pool
    with Pool(num_workers) as _pool:
        for _ in tqdm(
                _pool.imap_unordered(
                    __process_genome,
                    process_genome_arguments,
                ),
                ncols=80,
                smoothing=0.1,
                total=len(process_genome_arguments)):
            pass


# script (executable) for genome processing
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='process PATRIC genomes in parallel')

    parser.add_argument(
        '-i', '--genome_parent_dir_path', type=str,
        help='path to parent directory of all genomes)',
    )
    parser.add_argument(
        '-o', '--output_parent_dir_path', type=str, default=None,
        help='optional path to directory for processed genomes',
    )
    parser.add_argument(
        '-t', '--num_workers', type=int, default=1,
        help='number of workers for parallelization',
    )
    parser.add_argument(
        '-n', '--no_cd_search', action='store_true',
        help='disable the conserved domain searching')
    args = parser.parse_args()

    # usage example
    # $export PYTHONPATH=<project dir>:$PYTHONPATH
    # $python process_genomes.py \
    #      -i ~/data/PATRIC/genomes \
    #      -o ~/data/PATRIC/processed_genomes \
    #      -t 80
    process_genomes(
        genome_parent_dir_path=args.genome_parent_dir_path,
        output_parent_dir_path=args.output_parent_dir_path,
        num_workers=args.num_workers,
        no_cd_search=args.no_cd_search,
    )
