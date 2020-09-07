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
from Bio.Blast.Applications import NcbirpstblastnCommandline

from src import CDD_DIR_PATH

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


def conserved_domain_search(
        nucleotide_seq_path: str,
        cd_xml_path: Optional[str] = None,
) -> Blast:
    """conserved domain search for a given nucleotide sequence

    :param nucleotide_seq_path: path to nucleotide sequence to be searched
    for conserved domains
    :type nucleotide_seq_path: str
    :param cd_xml_path: optional path to the result in XML format,
    preferably ended with '.xml'
    :type cd_xml_path: str
    :return: Blast result from rpsblast/rpstblastn
    :rtype: Bio.Blast.Record.Blast
    """

    # return of the result XML already exists
    if os.path.exists(cd_xml_path):
        with open(cd_xml_path, 'r') as _fh:
            return NCBIXML.read(_fh)

    rpstblastn_cmd = NcbirpstblastnCommandline(
        query=nucleotide_seq_path,
        **RPSTBLASTN_KWARGS,
    )
    cd_xml, rpsblast_cmd_error_msg = rpstblastn_cmd()

    # write to result XML file if given
    if cd_xml_path:
        with open(cd_xml_path, 'w+') as _fh:
            _fh.write(cd_xml)

    return NCBIXML.read(io.StringIO(cd_xml))
