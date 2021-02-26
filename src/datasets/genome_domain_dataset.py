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
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Optional, Set, Tuple

import torch
import pandas as pd
from tqdm import tqdm
from torch.utils.data import Dataset

from src import RAW_DATA_DIR_PATH


BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH = os.path.join(
    RAW_DATA_DIR_PATH, 'genomes', 'bacteria.csv')
REF_N_REP_BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH = os.path.join(
    RAW_DATA_DIR_PATH, 'genomes', 'reference_or_representative_bacteria.csv')
_LOGGER = logging.getLogger(__name__)
resource.setrlimit(
    resource.RLIMIT_NOFILE,
    (4096, resource.getrlimit(resource.RLIMIT_NOFILE)[1]),
)

# Default token for contig sequences
DEFAULT_CODING_REGION_SEP_TOKEN: str = '<SEP>'
DEFAULT_CONTIG_BEGIN_TOKEN: str = '<BOCTG>'
DEFAULT_CONTIG_END_TOKEN: str = '<EOCTG>'


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
class ContigWithConservedDomains:
    """Data class for a genome contig, annotated with features (PATRIC or
    RefSeq) and the corresponding conserved domains.
    """
    genome_id: str
    genome_name: Optional[str]
    organism: Optional[Organism]
    ncbi_taxon_id: str
    annotation: Annotation
    contig_accession: str
    contig_feature_df: pd.DataFrame
    contig_conserved_domain_df: pd.DataFrame
    contig_feature_csv_file_path: str
    contig_conserved_domain_csv_file_path: str


def __get_organism_from_genome_name(genome_name: str) -> Organism:
    """Get the organism from the name of the genome by naive string parse,
    that is, if its name contains any of the organism strings, the genome
    is of that particular organism.

    :param genome_name:
    :type genome_name:
    :return:
    :rtype:
    """
    __lower_genome_name = genome_name.lower()
    for __organism in Organism:
        if __organism.value in __lower_genome_name:
            return __organism
    return Organism.OTHERS


def _get_contigs_with_conserved_domains(
        annotation: Annotation,
        genome_parent_dir_path: str,
        genome_summary_csv_file_path: Optional[str] = None,
) -> List[ContigWithConservedDomains]:
    """Get all the contigs inside a parent directory path into a list of
    ContigWithConservedDomains, which is essentially a data class
    with all the information on features and conserved domain annotations.

    Usage:
    contigs = _get_contigs_with_conserved_domains(
        annotation=Annotation.PATRIC,
        genome_parent_dir_path=
            '/home/xduan7/projects/bioseq-learning/data'
            '/interim/genomes/reference_or_representative_bacteria',
        genome_summary_csv_file_path=
            '/home/xduan7/projects/bioseq-learning'
            '/data/raw/genomes/reference_or_representative_bacteria.csv',
    )
    import pickle
    with open('reference_or_representative_bacteria_contigs_'
              'with_conserved_domains.PATRIC.pickle', 'wb') as _fh:
        pickle.dump(contigs, _fh)

    :param annotation:
    :type annotation:
    :param genome_parent_dir_path:
    :type genome_parent_dir_path:
    :param genome_summary_csv_file_path: optional genome summary CSV file
    path, which could be downloaded from PATRIC server. If given, this
    function will only process the genomes included in the CSV file by
    checking the "genome_id" column.
    :type genome_summary_csv_file_path: str
    :return:
    :rtype:
    """

    # load the genome summary dataframe
    try:
        _genome_summary_df = pd.read_csv(
            genome_summary_csv_file_path,
            index_col=None,
            usecols=[
                'Genome ID',
                'Genome Name',
                'Organism Name',
                'NCBI Taxon ID',
            ],
            dtype={
                'Genome ID': str,
                'Genome Name': str,
                'Organism Name': str,
                'NCBI Taxon ID': int,
            }
        )
        _genome_summary_df.columns = [
            'genome_id',
            'genome_name',
            'organism_name',
            'ncbi_taxon_id',
        ]
        _genome_summary_df = _genome_summary_df.set_index('genome_id')
        _genome_ids: Set[str] = set(_genome_summary_df.index.values)
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

    # construct conserved domain data class for every contig
    _contig_conserved_domains = []
    for __contig_conserved_domain_csv_file_path in \
            tqdm(_contig_conserved_domain_csv_file_paths):

        __split_path = \
            __contig_conserved_domain_csv_file_path.split(os.sep)
        __genome_id = __split_path[-3]
        __contig_feature_csv_file_path = \
            __contig_conserved_domain_csv_file_path.replace(
                '/conserved_domains/', '/features/').replace('.csv', '.tsv')

        # skip the config the the feature does not exist (should not happen)
        if not os.path.exists(__contig_feature_csv_file_path):
            _warning_msg = \
                f'The feature table file ({__contig_feature_csv_file_path}) ' \
                f'for current contig is missing. Skipping ...'
            _LOGGER.warning(_warning_msg)

        # skip the contig if the genome ID is not in the summary
        __genome_name, __organism = None, None
        if (_genome_summary_df is not None) and \
                (__genome_id not in _genome_ids):
            _warning_msg = \
                f'Genome {__genome_id} is not listed in the genome table ' \
                f'located in {genome_summary_csv_file_path}. Skipping ...'
            _LOGGER.warning(_warning_msg)
            continue

        __contig_accession = __split_path[-1].split('.', 1)[0]
        __contig_feature_df = pd.read_csv(
            __contig_feature_csv_file_path,
            sep='\t',
            header=0,
            index_col=None,
            dtype={'genome_id': str},
        )

        __contig_conserved_domain_df = pd.read_csv(
            __contig_conserved_domain_csv_file_path,
            header=0,
            index_col=None,
            dtype={
                'genome_id':            str,
                'pssm_id':              str,
                'superfamily_pssm_id':  str,
            }
        )

        # get the genome organism from genome name in the feature dataframe
        if __genome_name is None:
            __genome_names = __contig_feature_df.genome_name.unique().tolist()
            if len(__genome_names) > 1:
                __genome_name = max(__genome_names, key=len)
                _warning_msg = \
                    f'More than one genome names ({__genome_names}) in ' \
                    f'a single contig feature dataframe for contig ' \
                    f'{__contig_accession} in genome with ID {__genome_id}. ' \
                    f'Using the longest genome name {__genome_name} ...'
                _LOGGER.warning(_warning_msg)
            else:
                __genome_name = __genome_names[0]
            __organism = __get_organism_from_genome_name(__genome_name)

        # clean up the feature dataframe
        __contig_feature_df = __contig_feature_df[
            __contig_feature_df.accession == __contig_accession]
        if len(__contig_feature_df) == 0:
            _warning_msg = \
                f'There are no features for accession {__contig_accession} ' \
                f'in the feature table of genome with ID {__genome_id}.'
            _LOGGER.warning(_warning_msg)
            continue
        __contig_feature_df.drop('genome_id', axis=1, inplace=True)
        __contig_feature_df.drop('genome_name', axis=1, inplace=True)
        __contig_feature_df.drop('accession', axis=1, inplace=True)
        __contig_feature_df.drop('annotation', axis=1, inplace=True)

        # clean up the conserved domain dataframe
        __contig_conserved_domain_df.drop('genome_id', axis=1, inplace=True)
        __contig_conserved_domain_df.drop('genome_name', axis=1, inplace=True)

        __contig_with_conserved_domain = ContigWithConservedDomains(
            __genome_id,
            __genome_name,
            __organism,
            _genome_summary_df.loc[__genome_id, 'ncbi_taxon_id'],
            annotation,
            __contig_accession,
            __contig_feature_df,
            __contig_conserved_domain_df,
            __contig_feature_csv_file_path,
            __contig_conserved_domain_csv_file_path
        )
        _contig_conserved_domains.append(__contig_with_conserved_domain)
    return _contig_conserved_domains


def __convert_single_contig_with_conserved_domains_to_sequence(
        contig_with_conserved_domain: ContigWithConservedDomains,
        coding_region_sep_token: str = DEFAULT_CODING_REGION_SEP_TOKEN,
        contig_begin_token: Optional[str] = DEFAULT_CONTIG_BEGIN_TOKEN,
        contig_end_token: Optional[str] = DEFAULT_CONTIG_END_TOKEN,
) -> Tuple[str, List[str]]:

    _genome_id: str = contig_with_conserved_domain.genome_id






def _convert_contigs_with_conserved_domains_to_sequences(
        contig_conserved_domains: List[ContigWithConservedDomains],
        coding_region_sep_token: str = DEFAULT_CODING_REGION_SEP_TOKEN,
        contig_begin_token: Optional[str] = DEFAULT_CONTIG_BEGIN_TOKEN,
        contig_end_token: Optional[str] = DEFAULT_CONTIG_END_TOKEN,
) -> Dict[str, List[str]]:


    # get a list of arguments for
    # __convert_single_contig_with_conserved_domains_to_sequence
    __arg_list: List[Tuple[
        ContigWithConservedDomains,
        str,
        Optional[str],
        Optional[str],
    ]] = []





    # with Pool(cpu_count()) as _pool:
    #     _processed_contigs: \
    #         List[Optional[Tuple[str, int, np.ndarray, np.ndarray]]] = \
    #         _pool.starmap(
    #             __convert_single_contig_with_conserved_domains_to_sequence,
    #             __arg_list,
    #         )
    #
    #
    #
    # pass


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
