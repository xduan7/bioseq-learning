"""
File Name:          conserved_domain_search.py
Project:            bioseq-learning
File Description:
"""
import os
import re
import logging
from enum import Enum
from io import StringIO
from subprocess import Popen, PIPE
from typing import Any, Sequence, List, Optional, Union

import numpy as np
import pandas as pd
from Bio.Blast.Applications import \
    NcbirpsblastCommandline, \
    NcbirpstblastnCommandline, \
    NcbiblastformatterCommandline
from Bio.Application import ApplicationError

from src import CDD_DIR_PATH, CDD_DATA_DIR_PATH


_LOGGER = logging.getLogger(__name__)

# the following rpsblast/rpstblastn parameters are designed to replicate the
# CDD search results on NCBI website, using the entire CDD database and the
# BLAST archive (ASN.1) output format for direct processing with Biopython
# ref: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBFtp
RPSBLAST_KWARGS = {
    'db': CDD_DIR_PATH,
    'seg': 'no',
    'outfmt': 11,
    'evalue': 0.01,
    'comp_based_stats': '1',
}
# assuming rpsblast (for amino acid seq) and rpstblastn (for nucleotide seq)
# share the exact same arguments for conserved domain search
RPSTBLASTN_KWARGS = RPSBLAST_KWARGS

RPSBPROC_QUERY_COLUMNS = [
    '',
    'query_id',
    'seq_type',
    'seq_length',
    'fasta_comment',
]

# note that the domains columns are the same as superfamily's
RPSBPROC_DOMAIN_COLUMNS = [
    'session',
    'query_id[reading_frame]',
    'hit_type',
    'pssm_id',
    'start',
    'end',
    'e_value',
    'bitscore',
    'accession',
    'short_name',
    'incomplete',
    'superfamily_pssm_id',
]

PARSED_RPSBPROC_DOMAIN_COLUMNS = [
    'genome_id',
    'genome_name',
    'seq_id',
    'seq_description',
    'hit_type',
    'reading_frame',
    'short_name',
    'accession',
    'pssm_id',
    'start',
    'end',
    'e_value',
    'bitscore',
    'incomplete',
    'superfamily_short_name',
    'superfamily_accession',
    'superfamily_pssm_id',
]


class FASTAType(Enum):
    """type enum class for FASTA sequences
    *.fna: FASTA nucleic acid sequence for the whole genome
    *.faa: FASTA amino acid sequence of gene regions
    detailed reference: https://en.wikipedia.org/wiki/FASTA_format
    """
    fna = '.fna'
    faa = '.faa'


def __insert_column_value(
        target_dataframe: pd.DataFrame,
        target_column_name: str,
        source_seq: Sequence,
        source_seq_index: int,
        default_column_value: Any = '-',
):
    try:
        target_dataframe[target_column_name] = source_seq[source_seq_index]
    except IndexError:
        _info_msg = \
            f'cannot access index {source_seq_index} from source ' \
            f'sequence {source_seq} for column {target_column_name}.'
        _LOGGER.info(_info_msg)
        target_dataframe[target_column_name] = default_column_value
    except ValueError as __error:
        _info_msg = \
            f'cannot assign value {source_seq[source_seq_index]} to ' \
            f'column {target_column_name} due to \'{__error}\'.'
        _LOGGER.info(_info_msg)
        target_dataframe[target_column_name] = default_column_value


def __parse_rpsbproc_output(
        rpsbproc_output: str,
) -> pd.DataFrame:
    """helper function to parse the rpsbproc (post-processing for rpsblast)
    text, and return the result as a Pandas DataFrame
    :param rpsbproc_output:
    :type rpsbproc_output:
    :return:
    :rtype:
    """

    # remove all the comments in rpsbproc (starts with #)
    _rpsbproc_output = \
        re.sub(r'^#.*\n?', '', rpsbproc_output, flags=re.MULTILINE)

    # get the queries, domains, and the superfamilies
    # query is merely the header line for rpsbproc output
    _queries: List[str] = re.findall(
        r'\nQUERY(.*?)\nDOMAINS',
        _rpsbproc_output,
        re.DOTALL,
    )
    _domains: List[str] = re.findall(
        r'\nDOMAINS(.*?)\nENDDOMAINS',
        _rpsbproc_output,
        re.DOTALL,
    )
    _superfamilies: List[str] = re.findall(
        r'\nSUPERFAMILIES(.*?)\nENDSUPERFAMILIES',
        _rpsbproc_output,
        re.DOTALL,
    )
    # the number of queries should equal to the number of domains and
    # superfamilies; empty domains/superfamilies will be ''
    assert len(_queries) == len(_domains) == len(_superfamilies)

    _domain_superfamily_df = \
        pd.DataFrame(columns=PARSED_RPSBPROC_DOMAIN_COLUMNS)
    for __query, __domain, __superfamily in \
            zip(_queries, _domains, _superfamilies):

        # take the last part of query for FASTA definition
        __query_df = pd.read_csv(
            StringIO(__query),
            sep='\t',
            header=None,
        )
        __fasta_comment: str = \
            __query_df[__query_df.columns[-1]].to_list()[0]

        __domain_superfamily_df = pd.read_csv(
            StringIO(__domain + __superfamily),
            sep='\t',
            header=None,
            names=RPSBPROC_DOMAIN_COLUMNS,
        )

        # add information from FASTA comment
        # format: sequence id, sequence description [genome name | genome ID]
        __fasta_comment: List[str] = re.split(r' {2,}', __fasta_comment)
        __insert_column_value(
            target_dataframe=__domain_superfamily_df,
            target_column_name='seq_id',
            source_seq=__fasta_comment,
            source_seq_index=0,
        )
        __insert_column_value(
            target_dataframe=__domain_superfamily_df,
            target_column_name='seq_description',
            source_seq=__fasta_comment,
            source_seq_index=1,
        )
        __genome_info = __fasta_comment[2].strip('[]').split(' | ')
        __insert_column_value(
            target_dataframe=__domain_superfamily_df,
            target_column_name='genome_name',
            source_seq=__genome_info,
            source_seq_index=0,
        )
        __insert_column_value(
            target_dataframe=__domain_superfamily_df,
            target_column_name='genome_id',
            source_seq=__genome_info,
            source_seq_index=1,
        )

        # add reading frame column
        # if there is no reading frame, fill the cell with NaN
        __domain_superfamily_df['reading_frame'] = \
            __domain_superfamily_df['query_id[reading_frame]'].map(
                lambda _s:
                re.search(r'\[(.*?)\]', _s, re.DOTALL).group(1)
                if (('[' in _s) and (']' in _s)) else np.nan
            )

        # fix the superfamily column
        # 45 is not a valid superfamily PSSM ID, but the ASCII of '-'
        # make sure that all PSSM IDs are strings not integers
        __domain_superfamily_df['pssm_id'] = \
            __domain_superfamily_df['pssm_id'].apply(str)
        __domain_superfamily_df['superfamily_pssm_id'] = \
            __domain_superfamily_df['superfamily_pssm_id'].map(
                lambda _s: np.nan if (_s == '45' or _s == 45) else str(_s)
            )

        # add superfamily accession short name columns
        superfamily_df = __domain_superfamily_df.loc[
            __domain_superfamily_df['hit_type'] == 'Superfamily']
        superfamily_df = superfamily_df[
            ['short_name', 'accession', 'pssm_id']]
        superfamily_df.columns = [
            'superfamily_short_name',
            'superfamily_accession',
            'superfamily_pssm_id',
        ]
        __domain_superfamily_df = __domain_superfamily_df.merge(
            superfamily_df, 'outer', on='superfamily_pssm_id',
        ).fillna(np.nan)

        # merge the dataframe for the sequence
        # replace '-' for np.nan for empty/non-applicable cells
        _domain_superfamily_df = pd.concat([
            _domain_superfamily_df,
            __domain_superfamily_df.replace('-', np.nan),
        ])

    return _domain_superfamily_df[PARSED_RPSBPROC_DOMAIN_COLUMNS]


def search_conserved_domains(
        fasta_file_path: str,
        cd_ans_path: str,
        fasta_file_type: Optional[Union[FASTAType, str]] = None,
        cd_xml_path: Optional[str] = None,
        cd_txt_path: Optional[str] = None,
        cd_csv_path: Optional[str] = None,
) -> None:
    """perform conserved domain search for a given FASTA sequence and parse
    the results into multiple formats
    :param fasta_file_path: path to the FASTA file for domain search
    :type fasta_file_path: str
    :param fasta_file_type: FASTA file type, could be string or FASTAType
    defined in this file, or None, in which case the function will infer
    the FASTA type from the file extension
    :type fasta_file_type: Optional[Union[FASTAType, str]]
    :param cd_ans_path: path to the result in BLAST archive (ASN.1) format,
    preferably ended with '.ans'
    :type cd_ans_path: str
    :param cd_xml_path: optional path to the result in BLAST XML format,
    preferably ended with '.xml'
    :type cd_xml_path: Optional[str]
    :param cd_txt_path: optional path to the result in post-rpsblast in text
    format, preferably ended with '.txt'
    :type cd_txt_path: Optional[str]
    :param cd_csv_path: optional path to the result in post-rpsblast in CSV
    format, preferably ended with '.csv'
    :type cd_csv_path: Optional[str]
    :return: None
    """

    # refer the FASTA file type if not explicitly given or given as string
    if not fasta_file_type:
        fasta_file_type: str = os.path.splitext(fasta_file_path)[1]
    if not isinstance(fasta_file_type, FASTAType):
        try:
            fasta_file_type = FASTAType(fasta_file_type)
        except ValueError:
            _error_msg = \
                f'cannot parse the given FASTA file with extension ' \
                f'\'{fasta_file_type}\', which must be one of ' \
                f'{list(FASTAType.__members__.keys())}.'
            raise ValueError(_error_msg)

    # return of the result BLAST archive (ASN.1) already exists
    if not os.path.exists(cd_ans_path):
        if fasta_file_type == FASTAType.fna:
            rpsblast_cmd = NcbirpstblastnCommandline(
                query=fasta_file_path,
                **RPSTBLASTN_KWARGS,
            )
        elif fasta_file_type == FASTAType.faa:
            rpsblast_cmd = NcbirpsblastCommandline(
                query=fasta_file_path,
                **RPSTBLASTN_KWARGS,
            )
        else:
            _error_msg = \
                f'conserved domains search has not been implemented for ' \
                f'FASTA file type with extension \'{fasta_file_type.value}\'.'
            raise NotImplementedError(_error_msg)

        try:
            cd_ans, _ = rpsblast_cmd()
        except ApplicationError as __error:
            _warning_msg = f'error from rpsblast: {__error}; skipping ...'
            _LOGGER.warning(_warning_msg)
            return

        # write to result ANS.1 file if given
        if cd_ans_path:
            with open(cd_ans_path, 'w+') as _fh:
                _fh.write(cd_ans)

    # translate ANS to XML format for easier Biopython parsing
    if cd_xml_path and (not os.path.exists(cd_xml_path)):
        formatter_cmd = NcbiblastformatterCommandline(
            archive=cd_ans_path,
            out=cd_xml_path,
            outfmt=5,
        )
        _, formatter_cmd_error_msg = formatter_cmd()

    # post-rpsblast processing with rpsbproc and store in text format
    if cd_txt_path and (not os.path.exists(cd_txt_path)):
        rpsbproc_cmd = Popen(
            [
                f'rpsbproc',
                f'--infile', f'{cd_ans_path}',
                f'--outfile', f'{cd_txt_path}',
                f'--data-path', f'{CDD_DATA_DIR_PATH}',
                f'--data-mode', 'full',
                f'--evalue', f'{RPSTBLASTN_KWARGS["evalue"]}',
                f'--show-families',
                f'--quiet',
            ],
            stdout=PIPE,
            stderr=PIPE,
            stdin=PIPE,
        )
        rpsbproc_cmd.wait()
        rpsbproc_cmd.communicate()

    # parse the post-rpsblast processing results and store in CSV format
    if cd_csv_path and (not os.path.exists(cd_csv_path)):
        with open(cd_txt_path, 'r') as _fh:
            rpsbproc_output = _fh.read()
        rpsbproc_output_df = __parse_rpsbproc_output(rpsbproc_output)
        rpsbproc_output_df.to_csv(cd_csv_path, index=False)
