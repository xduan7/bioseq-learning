"""
File Name:          conserved_domain_search.py
Project:            bioseq-learning

File Description:

"""
import io
import os
from typing import Optional

from Bio.Blast import NCBIXML
from Bio.Blast.Record import Blast
from Bio.Blast.Applications import \
    NcbirpstblastnCommandline, NcbiblastformatterCommandline, \
    _NcbibaseblastCommandline, _Option

from src import CDD_DIR_PATH

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


def conserved_domain_search(
        nucleotide_seq_path: str,
        cd_ans_path: str,
        cd_xml_path: Optional[str],
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
