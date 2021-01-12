"""
File Name:          conserved_domain_search.py
Project:            bioseq-learning

File Description:

"""
import os
import re
from enum import Enum
from io import StringIO
from typing import Union, Optional
from subprocess import Popen, PIPE

import pandas as pd
from Bio.Blast.Applications import \
    NcbirpsblastCommandline, \
    NcbirpstblastnCommandline, \
    NcbiblastformatterCommandline

from src import CDD_DIR_PATH, CDD_DATA_DIR_PATH

os.environ['PATH'] += ':/home/xduan7/software/ncbi-blast/bin'

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


def parse_rpsbproc_output(
        rpsbproc_output: str,
) -> pd.DataFrame:

    # remove all the comments in rpsbproc (starts with #)
    rpsbproc_output = \
        re.sub(r'^#.*\n?', '', rpsbproc_output, flags=re.MULTILINE)

    # get the domains and superfamilies in the text
    try:
        domain_txt = re.search(
            'DOMAINS(.*?)ENDDOMAINS',
            rpsbproc_output,
            re.DOTALL,
        ).group(1)
    except AttributeError:
        domain_txt = ''

    try:
        superfamily_txt = re.search(
            'SUPERFAMILIES(.*?)ENDSUPERFAMILIES',
            rpsbproc_output,
            re.DOTALL,
        ).group(1)
    except AttributeError:
        superfamily_txt = ''

    try:
        domain_superfamily_df = pd.read_csv(
            StringIO(domain_txt + superfamily_txt),
            sep='\t',
            header=None,
            names=RPSBPROC_DOMAIN_COLUMNS,
        )

        # add reading frame column
        domain_superfamily_df['reading_frame'] = \
            domain_superfamily_df['query_id[reading_frame]'].map(
                lambda _s:
                re.search(r'\[(.*?)\]', _s, re.DOTALL).group(1)
                if (('[' in _s) and (']' in _s)) else '-'
            )

        # fix the superfamily column
        # 45 is not a valid superfamily PSSM ID, but the ASCII of '-'
        # make sure that all PSSM IDs are strings not integers
        domain_superfamily_df['pssm_id'] = \
            domain_superfamily_df['pssm_id'].apply(str)
        domain_superfamily_df['superfamily_pssm_id'] = \
            domain_superfamily_df['superfamily_pssm_id'].map(
                lambda _s: '-' if _s == '45' or _s == 45 else str(_s)
            )

        # add superfamily accession short name columns
        superfamily_df = domain_superfamily_df.loc[
            domain_superfamily_df['hit_type'] == 'Superfamily']
        superfamily_df = superfamily_df[
            ['short_name', 'accession', 'pssm_id']]
        superfamily_df.columns = [
            'superfamily_short_name',
            'superfamily_accession',
            'superfamily_pssm_id',
        ]
        domain_superfamily_df = domain_superfamily_df.merge(
            superfamily_df, 'outer', on='superfamily_pssm_id',
        ).fillna('-')

    except pd.errors.EmptyDataError:
        domain_superfamily_df = pd.DataFrame(
            columns=PARSED_RPSBPROC_DOMAIN_COLUMNS)

    return domain_superfamily_df[PARSED_RPSBPROC_DOMAIN_COLUMNS]


def search_conserved_domains(
        fasta_file_path: str,
        fasta_file_type: Optional[Union[FASTAType, str]],
        cd_ans_path: str,
        cd_xml_path: Optional[str],
        cd_txt_path: Optional[str],
        cd_csv_path: Optional[str],
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
        cd_ans, _ = rpsblast_cmd()

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
        rpsbproc_output_df = parse_rpsbproc_output(rpsbproc_output)
        rpsbproc_output_df.to_csv(cd_csv_path, index=False)


# temporary test code
search_conserved_domains(
    fasta_file_path='./9.55.PATRIC.faa',
    fasta_file_type=None,
    cd_ans_path='./9.55.PATRIC.faa.ans',
    cd_xml_path='./9.55.PATRIC.faa.xml',
    cd_txt_path='./9.55.PATRIC.faa.txt',
    cd_csv_path='./9.55.PATRIC.faa.csv',
)
