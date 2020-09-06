"""
File Name:          process_genomes.py
Project:            bioseq-learning-cd

File Description:

"""
import os
import json

from typing import Optional, Iterator

import pandas as pd
from Bio import SeqIO, SeqRecord

from src import CDD_DIR_PATH
from src.utilities import create_directory


# the following rpsblast/rpstblastn parameters are designed to replicate the
# CDD search results on NCBI website, using the entire CDD database and the
# XML output format for direct processing with Biopython
# ref: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBFtp
RPSTBLASTN_KWARGS = {
    'db': CDD_DIR_PATH,
    'seg': 'no',
    'outfmt': 5,
    'evalue': 0.01,
    'comp_based_stats': '1',
}


def process_genome(
        genome_dir_path: str,
        genome_id: Optional[str] = None,
):
    """process a PATRIC genome

    the processing steps are listed as the following:
    - prepare the directories
    - split all the contigs and their features
    - search for conserved domains on all the contigs
    - create an info file

    :param genome_dir_path: string path to PATRIC genome directory
    :param genome_id: optional string for genome ID; using the genome
    directory name if not given
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

    # create the subdirectories if not exist
    # - ./contigs: contig nucleotide sequences
    # - ./features:  feature annotations for every contig
    # - ./conserved_domains: conserved domain results for every contig
    contig_dir_path: str = os.path.join(genome_dir_path, 'contigs')
    feature_dir_path: str = os.path.join(genome_dir_path, 'features')
    conserved_domain_dir_path: str = \
        os.path.join(genome_dir_path, 'conserved_domains')
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
            os.path.join(feature_path, f'{_accession}.tsv')
        _source_mask = (_feature_by_accession_df['feature_type'] == 'source')
        _feature_by_accession_df[~_source_mask].to_csv(
            _feature_by_accession_path, index=None, sep='\t')

    # TODO: a helper function for conserved domain search and root tracing
    # TODO: conserved domain search for all contigs
    # TODO: include conserved domain information in contig_info_dict

    # save the genome information in json format
    genome_info_path: str = os.path.join(genome_dir_path, 'info.json')
    genome_info: dict = {
        'genome_name': feature_df['genome_name'].unique()[0],
        'number_of_contigs': len(contig_info_dict),
        'contigs': contig_info_dict,
    }
    with open(genome_info_path, 'w+') as _fh:
        json.dump(genome_info, _fh, indent=4)


def process_genomes(
        genome_parent_dir_path: str,
        num_threads: int = 1,
):
    """process PATRIC test_process_genomes in the given parent directory in parallel

    :param genome_parent_dir_path:
    :param num_threads:
    :return
    """
    pass
