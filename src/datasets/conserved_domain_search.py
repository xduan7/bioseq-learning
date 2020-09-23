"""
File Name:          conserved_domain_search.py
Project:            bioseq-learning

File Description:

"""
import os
import re
from io import StringIO
from typing import Optional

import pandas as pd
from Bio.Blast.Applications import \
    NcbirpstblastnCommandline, NcbiblastformatterCommandline

from src import CDD_DIR_PATH, CDD_DATA_DIR_PATH


# the following rpsblast/rpstblastn parameters are designed to replicate the
# CDD search results on NCBI website, using the entire CDD database and the
# BLAST archive (ASN.1) output format for direct processing with Biopython
# ref: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBFtp
RPSTBLASTN_KWARGS = {
    'db': CDD_DIR_PATH,
    'seg': 'no',
    'outfmt': 11,
    'evalue': 0.01,
    'comp_based_stats': '1',
}

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


def parse_rpsbproc_output(
        rpsbproc_output: str,
) -> pd.DataFrame:

    # remove all the comments in rpsbproc (starts with #)
    rpsbproc_output = re.sub(
        r'^#.*\n?', '', rpsbproc_output, flags=re.MULTILINE)

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


def conserved_domain_search(
        nucleotide_seq_path: str,
        cd_ans_path: str,
        cd_xml_path: Optional[str],
        cd_txt_path: Optional[str],
        cd_csv_path: Optional[str],
) -> None:
    """conserved domain search for a given nucleotide sequence

    :param nucleotide_seq_path: path to nucleotide sequence to be searched
    for conserved domains
    :type nucleotide_seq_path: str
    :param cd_ans_path: path to the result in BLAST archive (ASN.1) format,
    preferably ended with '.ans'
    :type cd_ans_path: str
    :param cd_xml_path: optional path to the result in BLAST XML format,
    preferably ended with '.xml'
    :type cd_xml_path: str
    :param cd_txt_path: optional path to the result in post-rpsblast in text
    format, preferably ended with '.txt'
    :type cd_txt_path: str
    :param cd_csv_path: optional path to the result in post-rpsblast in CSV
    format, preferably ended with '.csv'
    :type cd_csv_path: str
    :return: None
    """

    # return of the result BLAST archive (ASN.1) already exists
    if not os.path.exists(cd_ans_path):

        rpstblastn_cmd = NcbirpstblastnCommandline(
            query=nucleotide_seq_path,
            **RPSTBLASTN_KWARGS,
        )
        cd_ans, rpsblast_cmd_error_msg = rpstblastn_cmd()

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
        rpsbproc_cmd = os.popen(
            f'rpsbproc '
            f'--infile {cd_ans_path} '
            f'--outfile {cd_txt_path} '
            f'--data-path {CDD_DATA_DIR_PATH} '
            f'--data-mode full '
            f'--evalue {RPSTBLASTN_KWARGS["evalue"]} '
            # f'--show-families '
        )
        _ = rpsbproc_cmd.read()
        rpsbproc_cmd.close()

    # parse the post-rpsblast processing results and store in CSV format
    if cd_csv_path and (not os.path.exists(cd_csv_path)):
        with open(cd_txt_path, 'r') as _fh:
            rpsbproc_output = _fh.read()
        rpsbproc_output_df = parse_rpsbproc_output(rpsbproc_output)
        rpsbproc_output_df.to_csv(cd_csv_path, index=False)
