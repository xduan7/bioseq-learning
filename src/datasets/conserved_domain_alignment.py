"""
File Name:          conserved_domain_alignment.py
Project:            bioseq-learning

File Description:

"""
import os
import subprocess

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.decomposition import PCA

from src.utilities import create_directory
from src import INTERIM_DATA_DIR_PATH, PROCESSED_DATA_DIR_PATH



NUM_PCA_COMPONENT = 1024
ALGN_VALUE = ['seq_identity', 'bitscore', 'e_value'][0]


CDD_ID_PATH = os.path.join(INTERIM_DATA_DIR_PATH, 'CDD_metadata/cddid.tbl')

INTERIM_CDD_ALGN_DIR_PATH = \
    os.path.join(INTERIM_DATA_DIR_PATH, 'CDD_alignment')
create_directory(INTERIM_CDD_ALGN_DIR_PATH)

CDD_MASTER_SEQ_PATH = \
    os.path.join(INTERIM_CDD_ALGN_DIR_PATH, 'cdd_masters.fa')
CDD_MASTER_DB_PATH = \
    os.path.join(INTERIM_CDD_ALGN_DIR_PATH, 'cdd_masters.dmnd')
CDD_MASTER_INTERIM_ALGN_PATH = \
    os.path.join(INTERIM_CDD_ALGN_DIR_PATH, 'cdd_master_alignment.tsv')

PROCESSED_CDD_ALGN_DIR_PATH = \
    os.path.join(PROCESSED_DATA_DIR_PATH, 'CDD_alignment')
create_directory(PROCESSED_CDD_ALGN_DIR_PATH)

CDD_MASTER_PROCESSED_ALGN_PATH = \
    os.path.join(PROCESSED_CDD_ALGN_DIR_PATH, 'cdd_master_alignment.csv')
CDD_MASTER_PROCESSED_ALGN_FEAT_PATH = \
    os.path.join(PROCESSED_CDD_ALGN_DIR_PATH,
                 'cdd_master_alignment_feature.csv')


# read the identification metadata of CDs for PSSM <-> accession translation
cdd_id_df = pd.read_table(
    CDD_ID_PATH,
    sep='\t',
    header=None,
    names=[
        'pssm_id',
        'accession',
        'short_name',
        'description',
        'length'
    ],
    dtype={
        'pssm_id': str,
        'accession': str,
        'short_name': str,
        'description': str,
        'length': int,
    },
)

cdd_id_df = cdd_id_df.set_index('pssm_id')
pssm_id_to_accession_dict = cdd_id_df.to_dict()['accession']
cdd_master_seq_records = list(SeqIO.parse(CDD_MASTER_SEQ_PATH, 'fasta'))


# generate the databse and the pair-wise alignment of all CDs
if not os.path.exists(CDD_MASTER_DB_PATH):
    subprocess.run([
        'diamond', 'makedb',
        '--in', CDD_MASTER_SEQ_PATH,
        '--db', CDD_MASTER_DB_PATH,
        '--threads', '16',
    ])
if not os.path.exists(CDD_MASTER_INTERIM_ALGN_PATH):
    subprocess.run([
        'diamond', 'blastp',
        '--query', CDD_MASTER_SEQ_PATH,
        '--db', CDD_MASTER_DB_PATH,
        '--out', CDD_MASTER_INTERIM_ALGN_PATH,
        '--threads', '16',
        '--ultra-sensitive',
    ])

# rename and rearrange to columns of the CDD master sequence alignments
if not os.path.exists(CDD_MASTER_PROCESSED_ALGN_PATH):

    cdd_master_seq_algn_df = pd.read_csv(
        CDD_MASTER_INTERIM_ALGN_PATH,
        sep='\t',
        header=None,
        names=[
            'query_pssm_id',
            'target_pssm_id',
            'seq_identity',
            'length',
            'mismatches',
            'gap_openings',
            'query_start',
            'query_end',
            'target_start',
            'target_end',
            'e_value',
            'bitscore',
        ],
    )

    cdd_master_seq_algn_df['query_pssm_id'] = \
        cdd_master_seq_algn_df['query_pssm_id'].str.lstrip('gnl|CDD|')
    cdd_master_seq_algn_df['target_pssm_id'] = \
        cdd_master_seq_algn_df['target_pssm_id'].str.lstrip('gnl|CDD|')

    cdd_master_seq_algn_df['query_accession'] = \
        cdd_master_seq_algn_df['query_pssm_id'].map(pssm_id_to_accession_dict)
    cdd_master_seq_algn_df['target_accession'] = \
        cdd_master_seq_algn_df['target_pssm_id'].map(pssm_id_to_accession_dict)

    cdd_master_seq_algn_df = cdd_master_seq_algn_df[[
        'query_pssm_id', 'query_accession', 'target_pssm_id',
        'target_accession',
        'seq_identity', 'bitscore', 'e_value',
        'length', 'mismatches', 'gap_openings',
        'query_start', 'query_end', 'target_start', 'target_end',
    ]]

    cdd_master_seq_algn_df.to_csv(CDD_MASTER_PROCESSED_ALGN_PATH, index=False)
else:
    cdd_master_seq_algn_df = pd.read_csv(CDD_MASTER_PROCESSED_ALGN_PATH)

# # create symmetric matrix of sequence alignment for all CDs
# _cdd_master_seq_algn_df = cdd_master_seq_algn_df.loc[
#     cdd_master_seq_algn_df['query_pssm_id'] != cdd_master_seq_algn_df['target_pssm_id']]
# _cdd_master_seq_algn_df.columns = [
#     'target_pssm_id', 'target_accession', 'query_pssm_id', 'query_accession',
#     'seq_identity', 'bitscore', 'e_value',
#     'length', 'mismatches', 'gap_openings',
#     'target_start', 'target_end', 'query_start', 'query_end',
# ]
# _cdd_master_seq_algn_df = _cdd_master_seq_algn_df[[
#     'query_pssm_id', 'query_accession', 'target_pssm_id', 'target_accession',
#     'seq_identity', 'bitscore', 'e_value',
#     'length', 'mismatches', 'gap_openings',
#     'query_start', 'query_end', 'target_start', 'target_end',
# ]]

# # pivot the query and target accessions
# # this step causes int32 overflow error, need newer pandas version
# # reference: https://github.com/pandas-dev/pandas/issues/26314
# cdd_master_seq_algn_df = pd.concat(
#     [cdd_master_seq_algn_df, _cdd_master_seq_algn_df]
# ).pivot_table(
#     index='query_accession',
#     columns='target_accession',
#     values='seq_identity',
#     fill_value=0.,
#

accessions = list(set(
    list(cdd_master_seq_algn_df['query_accession'].unique()) +
    list(cdd_master_seq_algn_df['target_accession'].unique())
))

accession_index_dict = {
    _accession: _i for _i, _accession in enumerate(accessions)
}

cdd_master_seq_algn_mat = np.zeros(
    shape=(len(accession_index_dict), len(accession_index_dict)),
    dtype=np.float32,
)

for _row in cdd_master_seq_algn_df.itertuples():

    _query_index = accession_index_dict[_row.query_accession]
    _target_index = accession_index_dict[_row.target_accession]
    _algn_value = _row._asdict()[ALGN_VALUE]

    cdd_master_seq_algn_mat[_query_index, _target_index] = _algn_value
    cdd_master_seq_algn_mat[_target_index, _query_index] = _algn_value


# reduce the dimensionality and save for visualization
pca = PCA(n_components=NUM_PCA_COMPONENT)
pd.DataFrame(
    pca.fit_transform(cdd_master_seq_algn_mat),
    index=accessions,
).to_csv(CDD_MASTER_PROCESSED_ALGN_FEAT_PATH)
