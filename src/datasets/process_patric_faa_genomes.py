"""
File Name:          process_patric_faa_genomes.py
Project:            bioseq-learning

File Description:

"""
import os
import sys
import json
import logging
import argparse
import traceback
from multiprocessing import Pool
from typing import Optional, Iterator, List, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


_LOGGER = logging.getLogger(__name__)


def process_patric_faa_genome(
        genome_dir_path: str,
        output_dir_path: Optional[str] = None,
        genome_id: Optional[str] = None,
        no_cd_search: bool = False,
):
    """process a PATRIC *.faa genome (amino acid sequence)

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

    # iterate through the annotations (PATRIC and RefSeq)
    for _annotation in ['PATRIC', 'RefSeq']:

        # make sure that both the amino acid sequence AND the feature
        # table exist for the given annotation method
        _seq_path: str = os.path.join(
            genome_dir_path,
            f'{genome_id}.{_annotation}.faa',
        )
        _feature_df_path: str = os.path.join(
            genome_dir_path,
            f'{genome_id}.{_annotation}.features.tab',
        )

        if os.path.exists(_seq_path) and os.path.exists(_feature_df_path):

            # create the subdirectories if not exist
            # - ./OUTPUT_DIR: output directory if given
            # - ./OUTPUT_DIR/contigs: contig nucleotide sequences
            # - ./OUTPUT_DIR/features:  feature annotations for every contig
            # - ./OUTPUT_DIR/conserved_domains: conserved domains for every
            #   contig
            _output_dir_path: str = output_dir_path \
                if output_dir_path else genome_dir_path
            _contig_dir_path: str = os.path.join(output_dir_path, 'contigs')
            _feature_dir_path: str = os.path.join(output_dir_path, 'features')
            _conserved_domain_dir_path: str = \
                os.path.join(output_dir_path, 'conserved_domains')

            os.makedirs(_output_dir_path, exist_ok=True)
            os.makedirs(_contig_dir_path, exist_ok=True)
            os.makedirs(_feature_dir_path, exist_ok=True)
            os.makedirs(_conserved_domain_dir_path, exist_ok=True)

            # load the amino acid sequences for the whole genome
            # note that they are not continuous unless split into contigs
            _seq_record_list: List[SeqRecord] = \
                list(SeqIO.parse(_seq_path, 'fasta'))

            # load the feature table and split based on accession (contig)
            contig_info_dict: dict = {}
            _feature_df: pd.DataFrame = pd.read_table(_feature_df_path)
            for _accession, _accession_feature_df \
                    in _feature_df.groupby('accession'):

                # paths for the current accession
                _accession_feature_df_path: str = \
                    os.path.join(
                        _feature_dir_path,
                        f'{_accession}.{_annotation}.tsv'
                    )
                _accession_seq_path: str = \
                    os.path.join(
                        _contig_dir_path,
                        f'{_accession}.{_annotation}.faa'
                    )

                # save the feature tables for each accession (contig)
                _source_mask: pd.Series = \
                    (_accession_feature_df['feature_type'] == 'source')
                _accession_feature_df[~_source_mask].to_csv(
                    _accession_feature_df_path, index=None, sep='\t')

                # get all the amino acid sequences in the accession (contig)
                _accession_seq_records: List[SeqRecord] = []
                for __seq_record in _seq_record_list:
                    if __seq_record.id.count('|') == 1:
                        __seq_id = __seq_record.id
                    elif __seq_record.id.count('|') >= 2 \
                            and _annotation == 'PATRIC':
                        __seq_id = \
                            __seq_record.id.rstrip('|').rsplit('|', 1)[0]
                    elif __seq_record.id.count('|') >= 2 \
                            and _annotation == 'RefSeq':
                        __seq_id = \
                            __seq_record.id.rstrip('|').rsplit('|', 2)[1]
                    else:
                        _warning_msg = \
                            f'cannot parse the PATRIC ID from FASTA ' \
                            f'sequence record with name {__seq_record.id}.'
                        _LOGGER.warning(_warning_msg)
                        continue

                    # seq Id could either be PATRIC ID or RefSeq locus tag
                    # otherwise it does not belong to the current accession
                    if _annotation == 'PATRIC' and __seq_id in \
                            _accession_feature_df['patric_id'].values:
                        _accession_seq_records.append(__seq_record)
                    elif _annotation == 'RefSeq' and __seq_id in \
                            _accession_feature_df['refseq_locus_tag'].values:
                        _accession_seq_records.append(__seq_record)
                    else:
                        # this debug message is not necessary considering
                        # that we are analyzing the genome by
                        # contig/accession, and the sequence IDs are
                        # extracted from ALL the proteins on genome (not
                        # contig), which means that it's completely normal
                        # that we cannot find a sequence with its ID on the
                        # feature table of the current contig
                        # _debug_msg = \
                        #     f'faa sequence {__seq_id} in genome ' \
                        #     f'{genome_id} is not included in ' \
                        #     f'{_annotation} feature table.'
                        # _LOGGER.debug(_debug_msg)
                        continue

                if len(_accession_seq_records) > 0:
                    with open(_accession_seq_path, 'w') as __fh:
                        SeqIO.write(_accession_seq_records, __fh, 'fasta')
                    contig_info_dict[_accession] = {
                        'number_of_amino_acid_seqs':
                            len(_accession_seq_records)
                    }

            # iterate through all the contigs again for conserved domain search
            if not no_cd_search:

                # append the project source path to PYTHONPATH
                _project_src_dir = os.path.abspath(
                    os.path.join(__file__, '../../'))
                sys.path.append(_project_src_dir)
                from src.datasets.search_conserved_domains import \
                    search_conserved_domains

                for _contig_id in contig_info_dict.keys():
                    _contig_seq_path: str = os.path.join(
                        _contig_dir_path,
                        f'{_contig_id}.{_annotation}.faa',
                    )
                    _contig_cd_ans_path: str = os.path.join(
                        _conserved_domain_dir_path,
                        f'{_contig_id}.{_annotation}.ans',
                    )
                    _contig_cd_xml_path: str = os.path.join(
                        _conserved_domain_dir_path,
                        f'{_contig_id}.{_annotation}.xml',
                    )
                    _contig_cd_txt_path: str = os.path.join(
                        _conserved_domain_dir_path,
                        f'{_contig_id}.{_annotation}.txt',
                    )
                    _contig_cd_csv_path: str = os.path.join(
                        _conserved_domain_dir_path,
                        f'{_contig_id}.{_annotation}.csv',
                    )
                    search_conserved_domains(
                        fasta_file_path=_contig_seq_path,
                        cd_ans_path=_contig_cd_ans_path,
                        cd_xml_path=_contig_cd_xml_path,
                        cd_txt_path=_contig_cd_txt_path,
                        cd_csv_path=_contig_cd_csv_path,
                    )

            # save the genome information in json format
            try:
                genome_info_path: str = os.path.join(
                    output_dir_path,
                    f'info.{_annotation}.json',
                )
                genome_info: dict = {
                    'genome_name': _feature_df['genome_name'].unique()[0],
                    'number_of_contigs': len(contig_info_dict),
                    'contig_info': contig_info_dict,
                }
                with open(genome_info_path, 'w+') as _fh:
                    json.dump(genome_info, _fh, indent=4)
            except IndexError:
                _warning_msg = \
                    f'indexing error for genome {genome_id} with ' \
                    f'{_annotation} annotation. Please check the feature ' \
                    f'table {_feature_df_path}.'
                _LOGGER.warning(_warning_msg)

        # else:
        #     _info_msg = \
        #         f'Genome {genome_id} does not have complete {_annotation} ' \
        #         f'annotations, and therefore cannot be processed with the ' \
        #         f'annotated features and amino acid sequences.'
        #     _LOGGER.info(_info_msg)


def __process_patric_faa_genome(arg):
    """wrapper function for process_genome that only takes one tuple argument

    :param arg: tuple argument for process_genome to be unpacked
    :type arg: Tuple[str, Optional[str], Optional[str]]
    :return: None
    """
    try:
        return process_patric_faa_genome(*arg)
    except Exception as _e:
        _traceback_str = ''.join(traceback.format_tb(_e.__traceback__))
        _warning_msg = \
            f'failed to process the faa genome sequences in directory ' \
            f'\'{arg[0]}\' due to {_e.__class__.__name__}; traceback ' \
            f'of the exception:\n{_traceback_str}'
        _LOGGER.warning(_warning_msg)


def process_patric_faa_genomes(
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
    :param no_cd_search: optional flag indicating no conserved domain
    search, default set to False
    :type no_cd_search: bool
    :param num_workers: maximum number of workers for parallelization
    :type num_workers: int
    :return: None
    """

    # clamp the number of workers between range [1, number of CPU cores]
    num_workers: int = max(1, min(num_workers, os.cpu_count()))

    # get all the paths to genomes, and output paths if possible
    process_patric_faa_genome_args: \
        List[Tuple[str, Optional[str], str, bool]] = []
    for _genome_id in os.listdir(genome_parent_dir_path):
        _genome_dir_path: str = \
            os.path.join(genome_parent_dir_path, _genome_id)
        if os.path.isdir(_genome_dir_path):
            _output_dir_path: Optional[str] = \
                os.path.join(output_parent_dir_path, _genome_id) \
                if output_parent_dir_path else None
            process_patric_faa_genome_args.append(
                (_genome_dir_path, _output_dir_path, _genome_id, no_cd_search)
            )

    # parallel genome processing with pool
    with Pool(num_workers) as _pool:
        for _ in tqdm(
                _pool.imap_unordered(
                    __process_patric_faa_genome,
                    process_patric_faa_genome_args,
                ),
                ncols=80,
                smoothing=0.1,
                total=len(process_patric_faa_genome_args)):
            pass


# script (executable) for genome processing
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='process PATRIC *.faa genomes in parallel')

    parser.add_argument(
        '-i', '--genome_parent_dir_path', type=str,
        help='path to parent directory of all genomes)',
    )
    parser.add_argument(
        '-o', '--output_parent_dir_path', type=str, default=None,
        help='optional path to directory for processed genomes',
    )
    parser.add_argument(
        '-w', '--num_workers', type=int, default=1,
        help='number of workers for parallelization',
    )
    parser.add_argument(
        '-n', '--no_cd_search', action='store_true',
        help='disable the conserved domain searching',
    )
    args = parser.parse_args()

    # usage example
    # $export PYTHONPATH=<project dir>:$PYTHONPATH
    # $python process_genomes.py \
    #      -i ~/data/PATRIC/genomes \
    #      -o ~/data/PATRIC/processed_genomes \
    #      -w 80
    process_patric_faa_genomes(
        genome_parent_dir_path=args.genome_parent_dir_path,
        output_parent_dir_path=args.output_parent_dir_path,
        num_workers=args.num_workers,
        no_cd_search=args.no_cd_search,
    )
