"""
File Name:          genome_domain_dataset.py
Project:            bioseq-learning

File Description:

This file contains functions and classes for the conserved domain dataset
for genomes. Each eligible genome will be transformed into a sentence
of words of conserved domains, along with some special words.
TODO: what would be the target of each "sentence"?

The workflow of this dataset class initialization is shown below:
1. get all the contig conserved domain csv files
2. perform train/test split of csv files
    - randomly
    - stratified on organism (with hold-out organism)
3. get the vocabulary of domains from training set
4. construct training, validation, and test sets
    (1) tokenize the genome contigs with special characters
    (2) summarize the number of sequences and prepare the indexing

"""
import os
import logging
import resource
from glob import glob
from enum import Enum
from dataclasses import dataclass
from typing import Dict, List, Optional

import torch
import pandas as pd
from torch.utils.data import Dataset

from src import RAW_DATA_DIR_PATH


BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH = os.path.join(
    RAW_DATA_DIR_PATH, 'genomes', 'bacteria.csv')
_LOGGER = logging.getLogger(__name__)
resource.setrlimit(
    resource.RLIMIT_NOFILE,
    (4096, resource.getrlimit(resource.RLIMIT_NOFILE)[1]),
)


class Annotation(Enum):
    """Enum class for genome annotation sources in PATRIC database (
    reference: https://docs.patricbrc.org/user_guides/organisms_taxon/
    genome_annotations.html)
    """
    PATRIC = 'PATRIC'
    RefSeq = 'RefSeq'


class Organism(Enum):
    """Enum class for bacterial organisms in PATRIC database (reference:
    https://en.wikipedia.org/wiki/PATRIC).
    """
    BACILLUS        = 'bacillus'
    BARTONELLA      = 'bartonella'
    BORRELIA        = 'borrelia'
    BRUCELLA        = 'brucella'
    BURKHOLDERIA    = 'burkholderia'
    CAMPYLOBACTER   = 'campylobacter'
    CHLAMYDOPHILA   = 'chlamydophila'
    CLOSTRIDIUM     = 'clostridium'
    COXIELLA        = 'coxiella'
    EHRLICHIA       = 'ehrlichia'
    ESCHERICHIA     = 'escherichia'
    FRANCISELLA     = 'francisella'
    HELICOBACTER    = 'helicobacter'
    LISTERIA        = 'listeria'
    MYCOBACTERIUM   = 'mycobacterium'
    RICKETTSIA      = 'rickettsia'
    SALMONELLA      = 'salmonella'
    SHIGELLA        = 'shigella'
    STAPHYLOCOCCUS  = 'staphylococcus'
    VIBRIO          = 'vibrio'
    YERSINIA        = 'yersinia'
    OTHERS          = 'others'


@dataclass
class ContigConservedDomains:
    """Data class for a genome contig, annotated with conserved domains
    """
    genome_id: str
    genome_name: Optional[str]
    organism: Optional[Organism]
    annotation: Annotation
    contig_accession: str
    contig_conserved_domain_df: pd.DataFrame
    contig_conserved_domain_csv_file_path: str


@dataclass
class ContigConservedDomains:
    pass


def __get_organism_from_genome_name(genome_name: str) -> Organism:
    __lower_genome_name = genome_name.lower()
    for __organism in Organism:
        if __organism.value in __lower_genome_name:
            return __organism
    return Organism.OTHERS


def __get_contig_conserved_domains(
        annotation: Annotation,
        genome_parent_dir_path: str,
        genome_summary_csv_file_path: Optional[str] =
        BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH,
) -> List[ContigConservedDomains]:

    # load the genome summary dataframe
    try:
        _genome_summary_df = pd.read_csv(
            genome_summary_csv_file_path,
            index_col=None,
            usecols=['Genome ID', 'Genome Name'],
            dtype={'Genome ID': str, 'Genome Name': str}
        )
        _genome_summary_df.columns = ['genome_id', 'genome_name']
    except (ValueError, FileNotFoundError):
        _warning_msg = \
            f'Failed to load the summary dataframe for all the ' \
            f'genomes in directory {genome_parent_dir_path}.'
        _LOGGER.warning(_warning_msg)
        _genome_summary_df = None

    # get all the paths to the *.{annotation}.csv files in parent dir
    genome_parent_dir_path = os.path.abspath(genome_parent_dir_path)
    _contig_conserved_domain_csv_file_path_pattern = os.path.join(
        genome_parent_dir_path, '**', f'*.{annotation.value}.csv')
    _contig_conserved_domain_csv_file_paths: List[str] = glob(
        _contig_conserved_domain_csv_file_path_pattern, recursive=True)
    print(_contig_conserved_domain_csv_file_paths)

    # construct conserved domain data class for every contig
    _contig_conserved_domains = []
    for __contig_conserved_domain_csv_file_path in \
            _contig_conserved_domain_csv_file_paths:

        __split_path = \
            __contig_conserved_domain_csv_file_path.split(os.sep)
        __genome_id = __split_path[-3]

        # skip the contig if the genome ID is not in the summary
        if _genome_summary_df is None:
            __genome_name, __organism = None, None
        elif __genome_id not in _genome_summary_df['genome_id'].values:
            _warning_msg = \
                f'Genome {__genome_id} is not listed in the genome table ' \
                f'located in {genome_summary_csv_file_path}. Skipping ...'
            _LOGGER.warning(_warning_msg)
            continue
        else:
            __genome_info = _genome_summary_df[
                _genome_summary_df['genome_id'] == __genome_id]
            if len(__genome_info) > 1:
                _warning_msg = \
                    f'Multiple {len(__genome_info)} records in the genome ' \
                    f'table located in {genome_summary_csv_file_path} for ' \
                    f'genome {__genome_id}. Skipping ...'
                _LOGGER.warning(_warning_msg)
                continue
            __genome_name = __genome_info['genome_name'].tolist()[0]
            __organism = __get_organism_from_genome_name(__genome_name)

        __contig_accession = __split_path[-1].split('.', 1)[0]
        __contig_conserved_domain_df = pd.read_csv(
            __contig_conserved_domain_csv_file_path,
            header=0,
            index_col=None,
            dtype={
                'genome_id': str,
                'pssm_id': str,
                'superfamily_pssm_id': str,
            }
        )
        __contig_conserved_domain = ContigConservedDomains(
            __genome_id,
            __genome_name,
            __organism,
            annotation,
            __contig_accession,
            __contig_conserved_domain_df,
            __contig_conserved_domain_csv_file_path
        )
        _contig_conserved_domains.append(__contig_conserved_domain)
    return _contig_conserved_domains


class GenomeDomainDataset(Dataset):
    """Dataset class for (conserved) domains on genomes.
    """
    def __init__(
            self,
            vocab: Dict[str, int],
            seq_len: int,
            max_num_paddings: int,
    ):
        pass

    def __len__(self) -> int:
        pass

    def __getitem__(
            self,
            index: int,
    ) -> torch.LongTensor:
        # TODO: target for the sequence of genome domains (e.g. protein family)
        pass

# class GenomeDomainIterDataset(IterableDataset, GenomeDataset)
