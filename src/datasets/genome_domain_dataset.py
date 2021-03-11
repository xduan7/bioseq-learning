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
import pickle
import logging
import resource
from glob import glob
from enum import Enum
from bisect import bisect
from collections import Counter
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Optional, Set, Tuple

import torch
import pandas as pd
from tqdm import tqdm
from torch.utils.data import Dataset

from src import RAW_DATA_DIR_PATH, PROCESSED_DATA_DIR_PATH


BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH = os.path.join(
    RAW_DATA_DIR_PATH, 'genomes', 'bacteria.csv')
REF_N_REP_BACTERIA_GENOME_SUMMARY_CSV_FILE_PATH = os.path.join(
    RAW_DATA_DIR_PATH, 'genomes', 'reference_or_representative_bacteria.csv')
REF_OR_REP_BACTERIA_CONTIGS_WITH_CDS_FILE_PATH = os.path.join(
    PROCESSED_DATA_DIR_PATH,
    'genomes/reference_or_representative_bacteria_contigs_'
    'with_conserved_domains.{annotation}.pickle'
)
REF_OR_REP_BACTERIA_CONTIG_CDS_SEQS_FILE_PATH = os.path.join(
    PROCESSED_DATA_DIR_PATH,
    'genomes/reference_or_representative_bacteria_contig_'
    'sequences.{annotation}.pickle'
)

_LOGGER = logging.getLogger(__name__)
resource.setrlimit(
    resource.RLIMIT_NOFILE,
    (4096, resource.getrlimit(resource.RLIMIT_NOFILE)[1]),
)

# special makers for domain sequences
CONTIG_BEGIN_MARKER: str = '<bos>'
CONTIG_END_MARKER: str = '<eos>'
GENE_BEGIN_MARKER: str = '<cds>'
UNKNOWN_MARKER: str = '<unk>'
PADDING_MARKER: str = '<pad>'

SPECIAL_MARKERS: Set[str] = {
    GENE_BEGIN_MARKER,
    CONTIG_BEGIN_MARKER,
    CONTIG_END_MARKER,
    UNKNOWN_MARKER,
    PADDING_MARKER,
}
SPECIAL_MARKER_TOKENIZER: Dict[str, int] = {
    GENE_BEGIN_MARKER:   1,
    CONTIG_BEGIN_MARKER: 2,
    CONTIG_END_MARKER:   3,
    UNKNOWN_MARKER:      4,
    PADDING_MARKER:      0,
}


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


def __convert_single_contig_to_domain_sequence(
        contig_with_cds: ContigWithConservedDomains,
        include_superfamily: bool = True,
) -> Tuple[str, List[str]]:

    _id: str = \
        f'{contig_with_cds.genome_id}/' \
        f'{contig_with_cds.contig_accession}'
    _annotation: Annotation = contig_with_cds.annotation

    _feature_df: pd.DataFrame = \
        contig_with_cds.contig_feature_df
    _feature_df: pd.DataFrame = _feature_df[
        _feature_df['feature_type'] == 'CDS'
    ]
    if _annotation == Annotation.PATRIC:
        _feature_df: pd.DataFrame = _feature_df[
            ['patric_id', 'product', 'plfam_id', 'pgfam_id']]
    else:
        _feature_df: pd.DataFrame = _feature_df[
            ['refseq_locus_tag', 'product', 'plfam_id', 'pgfam_id']]
    _feature_df: pd.DataFrame = _feature_df.reset_index(drop=True)
    _feature_df.columns = ['seq_id', 'product', 'plfam_id', 'pgfam_id']

    _conserved_domain_df: pd.DataFrame = \
        contig_with_cds.contig_conserved_domain_df

    _hit_types = {'Specific', 'Superfamily'} \
        if include_superfamily else {{'Specific'}}
    _conserved_domain_df: pd.DataFrame = _conserved_domain_df[
        _conserved_domain_df['hit_type'].isin(_hit_types)
    ]
    _conserved_domain_df: pd.DataFrame = \
        _conserved_domain_df[[
            'seq_id', 'accession', 'hit_type', 'pssm_id',
            'start', 'end', 'e_value', 'bitscore',
    ]]
    _conserved_domain_df: pd.DataFrame = \
        _conserved_domain_df.reset_index(drop=True)

    def __get_seq_id(__seq_id: str):
        if __seq_id.count('|') == 1:
            return __seq_id
        elif __seq_id.count('|') >= 2 and \
                _annotation == Annotation.PATRIC:
            return __seq_id.rstrip('|').rsplit('|', 1)[0]
        elif __seq_id.count('|') >= 2 and \
                _annotation == Annotation.RefSeq:
            return __seq_id.rstrip('|').rsplit('|', 2)[1]
        else:
            _warning_msg = \
                f'cannot parse the PATRIC ID from FASTA ' \
                f'sequence record with name {__seq_id}.'
            print(_warning_msg)
            return ''

    _conserved_domain_df['seq_id'] = \
        _conserved_domain_df['seq_id'].apply(__get_seq_id)
    
    _feature_df = _feature_df.set_index('seq_id')
    _cds_seq_ids = _feature_df.index.values
    _ret_seq: List[str] = [CONTIG_BEGIN_MARKER]
    for __cds_seq_id in _cds_seq_ids:
        __cds_conserved_domain_df = _conserved_domain_df[
            _conserved_domain_df['seq_id'] == __cds_seq_id
        ].copy()
        __cds_conserved_domain_df.sort_values(
            by=['e_value', 'bitscore'],
            ascending=[True, False],
            inplace=True,
        )
        __cds_proc_conserved_domain_df = \
            pd.DataFrame([], columns=__cds_conserved_domain_df.columns)

        while len(__cds_conserved_domain_df) > 0:
            __curr_conserved_domain = \
                __cds_conserved_domain_df.iloc[0]
            __curr_start = __curr_conserved_domain.start
            __curr_end = __curr_conserved_domain.end
            __cds_proc_conserved_domain_df = \
                __cds_proc_conserved_domain_df.append(
                    __curr_conserved_domain,
                    ignore_index=True,
                )
            __cds_conserved_domain_df.drop(
                __cds_conserved_domain_df[(
                    (__cds_conserved_domain_df.start < __curr_end) &
                    (__cds_conserved_domain_df.end > __curr_start)
                )].index,
                inplace=True,
            )
            __cds_conserved_domain_df.reset_index(drop=True, inplace=True)

        __cds_proc_conserved_domain_df.sort_values(
            by=['start', 'bitscore'],
            inplace=True,
        )

        __cds_product = _feature_df.loc[__cds_seq_id, 'product']
        __cds_plfam_id = _feature_df.loc[__cds_seq_id, 'plfam_id']
        __cds_pgfam_id = _feature_df.loc[__cds_seq_id, 'pgfam_id']
        _ret_seq.append(f'{GENE_BEGIN_MARKER}/{__cds_product}/{__cds_plfam_id}/{__cds_pgfam_id}')
        _ret_seq.extend(__cds_proc_conserved_domain_df['accession'].to_list())

    _ret_seq.append(CONTIG_END_MARKER)
    return _id, _ret_seq


def _convert_contigs_to_domain_sequences(
        contigs_with_cds: List[ContigWithConservedDomains],
        include_superfamily: bool = True,
) -> Dict[str, List[str]]:
    __arg_list_for_single_contig: List[
        Tuple[ContigWithConservedDomains, bool]
    ] = []
    print('Preparing the arguments for contig conversion ...')
    for __contig_with_cds in tqdm(contigs_with_cds):
        __arg_list_for_single_contig.append((
            __contig_with_cds,
            include_superfamily,
        ))
    print('Converting contigs into sequences of conserved domains ...')
    with Pool(cpu_count()) as _pool:
        _contig_cds_seq: List[Tuple[str, List[str]]] = \
            _pool.starmap(
                __convert_single_contig_to_domain_sequence,
                tqdm(__arg_list_for_single_contig),
            )
    return {__c[0]: __c[1] for __c in _contig_cds_seq}


class GenomeDomainDataset(Dataset):
    """Dataset class for (conserved) domains on genomes.
    """
    __slots__ = [
        # arguments
        'annot',                        # Annotation
        'sized_seq_len',                # int
        'max_num_paddings',             # int
        'num_domain_vocab',             # int
        'num_gene_target_vocab',        # int
        'include_gene_target',          # bool
        # public attr
        'domain_seqs',                  # Dict[str, List[str]]
        'domain_vocab',                 # Set[str]
        'domain_tokenizer',             # Dict[str, int]
        'gene_target_vocab',            # Optional[Set[str]]
        'gene_target_tokenizer',        # Optional[Dict[str, int]]
        'tokenized_domain_seqs',        # Dict[str, torch.LongTensor]
        # private attr
        '_len'                          # int
        '_seq_ids'                      # Sequence[str]
        '_accumulated_num_sized_seqs'   # Sequence[int]
    ]

    def __init__(
            self,
            annot: Annotation,
            sized_seq_len: int,
            max_num_paddings: int,
            num_domain_vocab: int,
            num_gene_target_vocab: int,
    ):
        self.annot: Annotation = annot
        self.sized_seq_len: int = sized_seq_len
        self.max_num_paddings: int = max_num_paddings
        self.num_domain_vocab: int = num_domain_vocab
        self.num_gene_target_vocab: int = num_gene_target_vocab
        self.include_gene_target: bool = (num_gene_target_vocab >= 0)
        self._check_arg_sanity()

        self.domain_seqs: Dict[str, List[str]] = self._get_domain_seqs()
        self.domain_vocab: Set[str] = \
            self._get_domain_vocab()
        self.domain_tokenizer: Dict[str, int] = \
            self._get_domain_tokenizer()
        self.gene_target_vocab: Optional[Set[str]] = \
            self._get_gene_target_vocab()
        self.gene_target_tokenizer: Optional[Dict[str, int]] = \
            self._get_gene_target_tokenizer()
        self.tokenized_domain_seqs: Dict[str, torch.LongTensor] = \
            self._get_tokenized_domain_seqs()
        self._prepare_indexing()

    def _check_arg_sanity(self):
        if self.max_num_paddings > self.sized_seq_len:
            _warning_msg = \
                f'The maximum number of padding is greater than the ' \
                f'window size of the sequences, which could yield ' \
                f'sequences fully made of paddings.'
            _LOGGER.warning(_warning_msg)

    def _get_domain_seqs(self) -> Dict[str, List[str]]:
        contig_cds_seq_file_path = \
            REF_OR_REP_BACTERIA_CONTIG_CDS_SEQS_FILE_PATH.format(
                annotation=self.annot.value
            )
        if os.path.exists(contig_cds_seq_file_path):
            with open(contig_cds_seq_file_path, 'rb') as __fh:
                contig_cds_seqs: Dict[str, List[str]] = pickle.load(__fh)
        else:
            contigs_with_cds_file_path: str = \
                REF_OR_REP_BACTERIA_CONTIGS_WITH_CDS_FILE_PATH.format(
                    annotation=self.annot.value)
            with open(contigs_with_cds_file_path, 'rb') as __fh:
                contigs_with_cds: List[ContigWithConservedDomains] = \
                    pickle.load(__fh)
            contig_cds_seqs = _convert_contigs_to_domain_sequences(
                contigs_with_cds=contigs_with_cds,
            )
            with open(contig_cds_seq_file_path, 'wb') as __fh:
                pickle.dump(
                    contig_cds_seqs,
                    __fh, protocol=pickle.HIGHEST_PROTOCOL,
                )
        return contig_cds_seqs

    def _get_domain_vocab(self) -> Set[str]:
        if hasattr(self, 'domain_vocab'):
            return self.domain_vocab
        assert self.num_domain_vocab > len(SPECIAL_MARKERS)
        _domain_seqs: Dict[str, List[str]] = self._get_domain_seqs()
        _vocab_counter = Counter()
        for __seq in _domain_seqs.values():
            _vocab_counter.update(__seq)
        _vocab: Set[str] = SPECIAL_MARKERS.copy()
        for __v, _ in _vocab_counter.most_common():
            if __v.startswith(f'{GENE_BEGIN_MARKER}/'):
                continue
            _vocab.add(__v)
            if len(_vocab) == self.num_domain_vocab:
                break
        return _vocab

    def _get_domain_tokenizer(self) -> Dict[str, int]:
        if hasattr(self, 'domain_tokenizer'):
            return self.domain_tokenizer
        _token: Dict[str, int] = SPECIAL_MARKER_TOKENIZER
        for __v in self._get_domain_vocab():
            if __v not in _token:
                _token[__v] = len(_token)
        return _token

    def _get_gene_target_vocab(self) -> Optional[Set[str]]:
        if hasattr(self, 'gene_target_vocab'):
            return self.gene_target_vocab
        if not self.include_gene_target:
            return None
        _domain_seqs: Dict[str, List[str]] = self._get_domain_seqs()
        _vocab_counter = Counter()
        for __seq in _domain_seqs.values():
            _vocab_counter.update([
                # __d has the format of '<cds>/product/plfam_id/pgfam_id'
                # TODO: should we treat hypothetical proteins as if
                #  they are unknown?
                __d.split('/')[1] for __d in __seq
                if __d.startswith(f'{GENE_BEGIN_MARKER}/')
            ])
        print(f'total number of gene targets: {len(_vocab_counter)}')
        _vocab: Set[str] = set([UNKNOWN_MARKER] + [
            __v for __v, _ in
            _vocab_counter.most_common(self.num_gene_target_vocab - 1)
        ])
        return _vocab

    def _get_gene_target_tokenizer(self) -> Optional[Dict[str, int]]:
        if hasattr(self, 'gene_target_tokenizer'):
            return self.gene_target_tokenizer
        if not self.include_gene_target:
            return None
        _token: Dict[str, int] = {UNKNOWN_MARKER: 0}
        for __v in self._get_gene_target_vocab():
            if __v not in _token:
                _token[__v] = len(_token)
        return _token

    def _get_tokenized_domain_seqs(self) -> Dict[str, torch.LongTensor]:
        """get the tokenized domain sequences with padding
        """
        if hasattr(self, 'tokenized_cds_seqs'):
            return self.tokenized_cds_seqs
        print(
            f'Tokenizing {len(self.domain_seqs)} contig sequences '
            f'of domains with {len(self.domain_vocab)} tokens ...'
        )
        _unk_token: int = SPECIAL_MARKER_TOKENIZER[UNKNOWN_MARKER]
        _pad_token: int = SPECIAL_MARKER_TOKENIZER[PADDING_MARKER]
        _paddings: List[int] = [_pad_token] * self.max_num_paddings
        _tk_seqs: Dict[str, torch.LongTensor] = {}
        for __seq_id, __seq in tqdm(self._get_domain_seqs().items()):
            # split('/') to strip the gene targets from gene marker
            # e.g. '<cds>/product/plfam_id/pgfam_id' -> '<cds>'
            __tk_seq: torch.LongTensor = torch.LongTensor([
                self.domain_tokenizer.get(__cd.split('/')[0], _unk_token)
                for __cd in __seq
            ] + _paddings)
            _tk_seqs[__seq_id] = __tk_seq
        return _tk_seqs

    def _prepare_indexing(self):
        # sized sequences = padded sequences of given window size
        _seq_ids: List[str] = []
        _total_num_sized_domain_seqs: int = 0
        _accumulated_num_sized_domain_seqs: List[int] = []

        for __seq_id, __seq in self._get_tokenized_domain_seqs().items():
            _seq_ids.append(__seq_id)
            _num_seqs = len(__seq) - self.sized_seq_len + 1
            _total_num_sized_domain_seqs += _num_seqs
            _accumulated_num_sized_domain_seqs.append(
                _total_num_sized_domain_seqs)

        self._len = _total_num_sized_domain_seqs
        self._seq_ids = _seq_ids
        self._accumulated_num_sized_domain_seqs = \
            _accumulated_num_sized_domain_seqs

    def __len__(self) -> int:
        if not hasattr(self, '_len'):
            self._prepare_indexing()
        return self._len

    def __get_sized_tk_gene_target_seq(
            self,
            seq_id: str,
            sized_seq_start_pos: int,
            sized_seq_end_pos: int,
    ) -> Optional[torch.LongTensor]:
        if not self.include_gene_target:
            return None
        _domain_seq: List[str] = self.domain_seqs[seq_id]
        _sized_domain_seq: List[str] = \
            _domain_seq[sized_seq_start_pos: sized_seq_end_pos]
        _tk_gene_target_seq: List[int] = []
        _unk_token: int = self.gene_target_tokenizer[UNKNOWN_MARKER]
        for __d in _sized_domain_seq:
            if __d.startswith(f'{GENE_BEGIN_MARKER}/'):
                __gene_target: str = __d.split('/')[1]
                __tk_gene_target: int = \
                    self.gene_target_tokenizer.get(__gene_target, _unk_token)
                _tk_gene_target_seq.append(__tk_gene_target)
        return torch.LongTensor(_tk_gene_target_seq)

    def __getitem__(
            self,
            index: int,
    ) -> Tuple[torch.LongTensor, Optional[torch.LongTensor]]:
        # index = index of the data set
        # seq index = index of all contig sequences of domains
        _seq_index: int = \
            bisect(self._accumulated_num_sized_domain_seqs, index)
        _seq_id: str = self._seq_ids[_seq_index]

        _sized_seq_start_pos: int = index if _seq_index == 0 else \
            (index - self._accumulated_num_sized_domain_seqs[_seq_index - 1])
        _sized_seq_end_pos: int = _sized_seq_start_pos + self.sized_seq_len

        _tk_domain_seq: torch.LongTensor = \
            self.tokenized_domain_seqs[_seq_id]
        _sized_tk_domain_seq: torch.LongTensor = \
            _tk_domain_seq[_sized_seq_start_pos: _sized_seq_end_pos]

        _sized_tk_gene_target_id_seq: Optional[torch.LongTensor] = \
            self.__get_sized_tk_gene_target_seq(
                seq_id=_seq_id,
                sized_seq_start_pos=_sized_seq_start_pos,
                sized_seq_end_pos=_sized_seq_end_pos,
            )
        return _sized_tk_domain_seq, _sized_tk_gene_target_id_seq


dset = GenomeDomainDataset(
    annot=Annotation.PATRIC,
    sized_seq_len=10,
    max_num_paddings=5,
    num_domain_vocab=50,
    num_gene_target_vocab=50,
)

# class GenomeDomainIterDataset(IterableDataset, GenomeDataset)
