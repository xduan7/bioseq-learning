"""
File Name:          genome_visualization.py
Project:            bioseq-learning

File Description:

"""
import os
from typing import Dict, List

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqRecord

from src import DOC_DIR_PATH, E_COLI_GENOME_PARENT_DIR_PATH, \
    INTERIM_DATA_DIR_PATH, PROCESSED_DATA_DIR_PATH

image_dir_path = os.path.join(os.path.join(DOC_DIR_PATH, 'images'))
genome_dir_paths: List[str] = [
    os.path.join(E_COLI_GENOME_PARENT_DIR_PATH, _genome_id)
    for _genome_id in os.listdir(E_COLI_GENOME_PARENT_DIR_PATH)
    if os.path.isdir(os.path.join(E_COLI_GENOME_PARENT_DIR_PATH, _genome_id))
]


# number of contig count plot (for each genome) and contig length histogram
contig_cnt_hist_path = os.path.join(image_dir_path, 'contig_cnt_hist.png')
contig_len_hist_path = os.path.join(image_dir_path, 'contig_len_hist.png')
if (not os.path.exists(contig_cnt_hist_path)) or \
        (not os.path.exists(contig_len_hist_path)):

    genome_contig_len_dict: Dict[str, Dict[str, int]] = {}
    for __genome_dir_path in genome_dir_paths:

        __genome_id: str = os.path.basename(__genome_dir_path.rstrip('/'))
        __genome_contig_seq_dir_path: str = \
            os.path.join(__genome_dir_path, 'contigs')

        for __contig_seq_file_name in os.listdir(__genome_contig_seq_dir_path):
            __contig_seq_path = os.path.join(
                __genome_contig_seq_dir_path,
                __contig_seq_file_name,
            )

            with open(__contig_seq_path, 'r') as _fh:
                __contig_seq_rec: SeqRecord = next(SeqIO.parse(_fh, 'fasta'))

            __contig_id: str = f'{__contig_seq_rec.id}'
            __contig_seq: str = str(__contig_seq_rec.seq).lower()

            # count the contig only if the nucleotide bases are ATGC
            if len(set(__contig_seq) - {'a', 't', 'g', 'c'}) == 0:
                if __genome_id not in genome_contig_len_dict:
                    genome_contig_len_dict[__genome_id] = {}
                genome_contig_len_dict[__genome_id][__contig_id] \
                    = len(__contig_seq)

    genome_num_contig_list = [
        [__genome_id, len(__genome_contig_dict)]
        for __genome_id, __genome_contig_dict in genome_contig_len_dict.items()
    ]
    genome_num_contig_df = pd.DataFrame(
        genome_num_contig_list,
        columns=['genome_id', 'num_contigs'],
    )
    __genome_num_contig_df = genome_num_contig_df['num_contigs'].value_counts()
    for __i in range(1, 21):
        if __i not in __genome_num_contig_df:
            __genome_num_contig_df[__i] = 0
    plt.figure(figsize=(16, 9))
    sns.barplot(
        x=__genome_num_contig_df.index,
        y=__genome_num_contig_df,
        color='steelblue',
        edgecolor=None,
        saturation=1,
    )
    plt.title(
        f'contig count plot for e. coli '
        f'{len(genome_num_contig_df)} genomes'
    )
    plt.xlabel('number of contigs for a genome')
    plt.ylabel('number of genomes')
    plt.savefig(contig_cnt_hist_path)

    contig_len_list = []
    for __genome_id, __genome_contig_dict in genome_contig_len_dict.items():
        for __contig_id, __contig_seq_len in __genome_contig_dict.items():
            contig_len_list.append(
                [f'__genome_id/__contig_id', __contig_seq_len]
            )
    contig_len_df = pd.DataFrame(
        contig_len_list,
        columns=['genome_contig_id', 'contig_seq_len'],
    )
    plt.figure(figsize=(16, 9))
    sns.histplot(
        data=contig_len_df,
        x='contig_seq_len',
        bins=60,
        color='steelblue',
        edgecolor=None,
        alpha=1,
    )
    plt.xlim(0, 6e6)
    plt.title(f'length histogram plot for e. coli {len(contig_len_df)} contigs')
    plt.xlabel('contig length')
    plt.ylabel('number of contigs')
    plt.savefig(contig_len_hist_path)


# coding region (gene) histogram in terms of length
gene_len_hist_path = os.path.join(image_dir_path, 'gene_len_hist.png')
if not os.path.exists(gene_len_hist_path):
    
    __cds_len_list = []
    for __genome_dir_path in genome_dir_paths:

        __genome_contig_features_dir_path: str = \
            os.path.join(__genome_dir_path, 'features')

        for __contig_feature_file_name in \
                os.listdir(__genome_contig_features_dir_path):
            __contig_features_path = os.path.join(
                __genome_contig_features_dir_path,
                __contig_feature_file_name,
            )
            __contig_feature_df = pd.read_csv(__contig_features_path, sep='\t')
            # __contig_feature_df = __contig_feature_df[
            #     __contig_feature_df['plfam_id'].notna()]
            # __contig_feature_df = __contig_feature_df[
            #     __contig_feature_df['pgfam_id'].notna()]
            # __contig_feature_df = __contig_feature_df.loc[
            #     __contig_feature_df['product'] != ' hypothetical protein']
            __cds_len_list.extend(__contig_feature_df['na_length'].to_list())

    __max_cds_len = 6000
    __original_num_cds = len(__cds_len_list)
    __cds_len_list = [
        __cds_len for __cds_len in __cds_len_list
        if __cds_len <= __max_cds_len
    ]
    plt.figure(figsize=(16, 9))
    sns.histplot(
        x=__cds_len_list,
        bins=100,
        color='steelblue',
        edgecolor=None,
        alpha=1,
    )
    plt.xlim(0, __max_cds_len)
    plt.title(
        f'length histogram plot for {len(__cds_len_list)} CDS '
        f'with less than {__max_cds_len} bases ('
        f'{100.*len(__cds_len_list) / __original_num_cds:.2f}%)'
    )
    plt.xlabel('CDS length')
    plt.ylabel('number of CDS')
    plt.savefig(gene_len_hist_path)

# non-coding region histogram in terms of length
non_coding_len_hist_path = \
    os.path.join(image_dir_path, 'non_coding_len_hist.png')
if not os.path.exists(non_coding_len_hist_path):
    __non_coding_len_list = []
    for __genome_dir_path in genome_dir_paths:
        __genome_contig_features_dir_path: str = \
            os.path.join(__genome_dir_path, 'features')

        for __contig_feature_file_name in \
                os.listdir(__genome_contig_features_dir_path):
            __contig_features_path = os.path.join(
                __genome_contig_features_dir_path,
                __contig_feature_file_name,
            )
            __contig_feature_df = pd.read_csv(__contig_features_path, sep='\t')

            __start = __contig_feature_df['start'].to_list()
            __end = __contig_feature_df['end'].to_list()

            __non_coding_len_list.extend([
                (__start[__i + 1] - __end[__i])
                for __i in range(len(__start) - 1)
                if (__start[__i + 1] - __end[__i]) > 0
            ])


    plt.figure(figsize=(16, 9))
    sns.histplot(
        x=__non_coding_len_list,
        bins=100,
        color='steelblue',
        edgecolor=None,
        alpha=1,
    )
    plt.xlim(0, max(__non_coding_len_list))
    plt.title(
        f'length histogram plot for '
        f'{len(__non_coding_len_list)} non-coding regions'
    )
    plt.xlabel('non-coding region length')
    plt.ylabel('number of non-coding regions')
    plt.savefig(non_coding_len_hist_path)


# create conserved domain hits
cd_hit_csv_path = os.path.join(
    PROCESSED_DATA_DIR_PATH, 'genomes', 'cd_hit.csv')
if not os.path.exists(cd_hit_csv_path):
    CDD_SUPERFAMILY_PATH = os.path.join(
        INTERIM_DATA_DIR_PATH, 'CDD_metadata/family_superfamily_links')
    cdd_superfamily_df = pd.read_table(
        CDD_SUPERFAMILY_PATH,
        header=None,
        names=[
            'accession',
            'pssm_id',
            'superfamily_accession',
            'superfamily_pssm_id'
        ],
        dtype={
            'accession': str,
            'pssm_id': str,
            'superfamily_accession': str,
            'superfamily_pssm_id': str,
        },
    )
    cdd_superfamily_df['superfamily'] = \
        (cdd_superfamily_df['pssm_id'] ==
         cdd_superfamily_df['superfamily_pssm_id'])

    CDD_ID_PATH = os.path.join(INTERIM_DATA_DIR_PATH, 'CDD_metadata/cddid.tbl')
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

    cd_hit_df = pd.merge(
        cdd_id_df,
        cdd_superfamily_df[['accession', 'superfamily']],
        how='outer',
        on=['accession'],
    )
    cd_hit_df = cd_hit_df.fillna(True)
    cd_hit_df = cd_hit_df.set_index('pssm_id')
    cd_hit_df['length'] = cd_hit_df['length'] * 3
    cd_hit_df[['specific_hits', 'nonspecific_hits', 'other_hits']] = 0

    os.makedirs(os.path.dirname(cd_hit_csv_path), exist_ok=True)
    __unkown_pssm_ids = set()

    for __genome_dir_path in genome_dir_paths:

        __genome_contig_cd_dir_path: str = \
            os.path.join(__genome_dir_path, 'conserved_domains')

        for __contig_cd_file_name in os.listdir(__genome_contig_cd_dir_path):
            if not __contig_cd_file_name.endswith('.csv'):
                continue

            __contig_cd_file_path = os.path.join(
                __genome_contig_cd_dir_path,
                __contig_cd_file_name,
            )
            __contig_cd_df = pd.read_csv(
                __contig_cd_file_path,
                usecols=['pssm_id', 'accession', 'hit_type'],
                dtype={
                    'pssm_id': str,
                    'accession': str,
                    'hit_type': str,
                },
            )
            for __hit in __contig_cd_df.itertuples():

                if __hit.pssm_id not in cd_hit_df.index:
                    __unkown_pssm_ids.add(__hit.pssm_id)
                    continue

                if __hit.hit_type == 'Specific':
                    cd_hit_df.at[__hit.pssm_id, 'specific_hits'] \
                        = 1 + cd_hit_df.at[__hit.pssm_id, 'specific_hits']
                elif __hit.hit_type == 'Non-specific':
                    cd_hit_df.at[__hit.pssm_id, 'nonspecific_hits'] \
                        = 1 + cd_hit_df.at[__hit.pssm_id, 'nonspecific_hits']
                else:
                    cd_hit_df.at[__hit.pssm_id, 'other_hits'] \
                        = 1 + cd_hit_df.at[__hit.pssm_id, 'other_hits']

    print(f'{len(__unkown_pssm_ids)} PSSM IDs not found: {__unkown_pssm_ids}')
    cd_hit_df.to_csv(cd_hit_csv_path)


# conserved domain histogram in terms of length
cd_len_hist_path = os.path.join(image_dir_path, 'cd_len_hist.png')
if not os.path.exists(cd_len_hist_path):
    cd_hit_df = pd.read_csv(
        cd_hit_csv_path,
        index_col='pssm_id',
        dtype={
            'pssm_id': str,
            'accession': str,
            'short_name': str,
            'description': str,
            'length': int,
            'superfamily': bool,
            'specific_hits': int,
            'nonspecific_hits': int,
            'other_hits': int,
        },
    )
    cd_hit_df.index = cd_hit_df.index.astype(str)
    cd_hit_df['total_hits'] = (
            cd_hit_df['specific_hits'] +
            cd_hit_df['nonspecific_hits'] +
            cd_hit_df['other_hits']
    )
    cd_hit_df['has_hits'] = (cd_hit_df['total_hits'] > 0)

    __hit_percentage = 100 * len(cd_hit_df[cd_hit_df["has_hits"]]) / len(cd_hit_df)
    print(f'{len(cd_hit_df[cd_hit_df["has_hits"]])} out of {len(cd_hit_df)} '
          f'({__hit_percentage:2.3f}%) conserved domains are found '
          f'in {len(genome_dir_paths)} high quality E coli genomes ...')
    __super_family_hit_percentage = \
        100 * len(cd_hit_df[cd_hit_df["has_hits"] & cd_hit_df["superfamily"]]) \
        / len(cd_hit_df[cd_hit_df["superfamily"]])
    print(f'{len(cd_hit_df[cd_hit_df["has_hits"] & cd_hit_df["superfamily"]])} '
          f'out of {len(cd_hit_df[cd_hit_df["superfamily"]])} ('
          f'{__super_family_hit_percentage:2.3f}%) conserved domains '
          f'superfamilies are found in {len(genome_dir_paths)} '
          f'high quality E coli genomes ...')

    plt.figure(figsize=(16, 9))
    _max_cd_len = 5000
    sns.histplot(
        cd_hit_df[cd_hit_df['length'] <= _max_cd_len],
        x='length',
        hue='has_hits',
        multiple='stack',
        bins=60,
        color='steelblue',
        edgecolor=None,
        alpha=1,
    )
    plt.legend(title='in e coli genomes')
    plt.xlim(0)
    __num_cds = len(cd_hit_df[cd_hit_df['length'] <= _max_cd_len])
    plt.title(
        f'length histogram plot for {__num_cds} conserved domains '
        f'with length less than {_max_cd_len} (out of {len(cd_hit_df)} CDs)'
    )
    plt.xlabel('conserved domain length')
    plt.ylabel('number of conserved domains')
    plt.savefig(cd_len_hist_path)


# conserved domain histogram in terms of occurrence
# conserved domain root histogram in terms of occurrences
