"""
File Name:          genome_parsing.py
Project:            bioseq-learning-cd

File Description:

"""
# imports, constants for paths and configurations

import os
import time

from Bio import SeqIO
from Bio.Blast.Applications import NcbirpsblastCommandline


# path to CDD (conserved domain database) and related info
CDD_PATH = '../data/interim/CDD/cdd'
RPSBLAST_KWARGS = {
    'db': CDD_PATH,
    'seg': 'no',
    'comp_based_stats': '1',
    'evalue': 0.01,
    'outfmt': 5,
}

# path to all the genomes (downloaded from PATRIC directly)
GENOME_PARENT_PATH = './genome_parsing_examples'


# get all the genome names (PATRIC ID)
genome_names = next(os.walk(GENOME_PARENT_PATH))[1]


# loop over all the genomes on parent path and perform parsing ...
for _genome_name in genome_names:

    _genome_path = os.path.join(GENOME_PARENT_PATH, _genome_name)

    # get the nucleotide sequences (whole contig and the one with features only)
    _cntg_path = os.path.join(_genome_path, _genome_name + '.fna')
    _cntg_records = list(SeqIO.parse(_cntg_path, 'fasta'))
    _feat_path = os.path.join(_genome_path, _genome_name + '.PATRIC.ffn')
    _feat_records = list(SeqIO.parse(_feat_path, 'fasta'))

    print(f'Genome {_genome_name} contains '
          f'{len(_cntg_records)} contig(s) and '
          f'{len(_feat_records)} features (genes, RNAs, etc.)')

    # perform CD search for the contigs with measured time
    _start_time = time.time()
    _cntg_rpsblast_cmd = NcbirpsblastCommandline(query=_cntg_path, **RPSBLAST_KWARGS)
    _cntg_rpsblast_xml_result, _cntg_rpsblast_cmd_error_msg = _cntg_rpsblast_cmd()
    print(time.time() - _start_time)
    with open(os.path.join(GENOME_PARENT_PATH, 'rpsblast_contig_' + _genome_name + '.xml'), 'w+') as f:
        f.write(_cntg_rpsblast_xml_result)

    # perform CD search for the features only with measure time
    _start_time = time.time()
    _feat_rpsblast_cmd = NcbirpsblastCommandline(query=_feat_path, **RPSBLAST_KWARGS)
    _feat_rpsblast_xml_result, _feat_rpsblast_cmd_error_msg = _feat_rpsblast_cmd()
    print(time.time() - _start_time)
    with open(os.path.join(GENOME_PARENT_PATH, 'rpsblast_feature_' + _genome_name + '.xml'), 'w+') as f:
        f.write(_feat_rpsblast_xml_result)


# TODO: make sure that the two CD searches have the same results ...
