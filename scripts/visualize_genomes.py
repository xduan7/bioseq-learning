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

from src import DOC_DIR_PATH, E_COLI_GENOME_PARENT_DIR_PATH

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
gene_len_hist_path = os.path.join(image_dir_path, 'gene_cnt_hist.png')
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


# conserved domain histogram in terms of length
# conserved domain histogram in terms of occurrence
# conserved domain root histogram in terms of occurrences
