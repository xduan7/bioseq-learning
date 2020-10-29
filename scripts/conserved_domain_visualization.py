"""
File Name:          conserved_domain_visualization.py
Project:            bioseq-learning

File Description:

"""
import os
import logging

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from src import INTERIM_DATA_DIR_PATH, PROCESSED_DATA_DIR_PATH, DOC_DIR_PATH


MIN_SUPERFAMILY_SIZE = 100

CDD_SUPERFAMILY_PATH = os.path.join(
    INTERIM_DATA_DIR_PATH, 'CDD_metadata/family_superfamily_links')
CDD_MASTER_PROCESSED_PAIRWISE_ALGN_PATH = os.path.join(
    PROCESSED_DATA_DIR_PATH, 'CDD_alignment/cdd_master_pairwise_alignment.npy')

CDD_TSNE_IMAGE_DIR_PATH = os.path.join(DOC_DIR_PATH, 'images/cdd_tsne')
CDD_UMAP_IMAGE_DIR_PATH = os.path.join(DOC_DIR_PATH, 'images/cdd_umap')


numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)
matplotlib_logger = logging.getLogger('matplotlib')
matplotlib_logger.setLevel(logging.WARNING)

os.makedirs(CDD_TSNE_IMAGE_DIR_PATH, exist_ok=True)
os.makedirs(CDD_UMAP_IMAGE_DIR_PATH, exist_ok=True)

cdd_superfamily_df = pd.read_table(
    CDD_SUPERFAMILY_PATH,
    header=None,
    index_col=0,
    names=[
        'accession',
        'pssm_id',
        'superfamily_accession',
        'superfamily_pssm_id'
    ],
)

with open(CDD_MASTER_PROCESSED_PAIRWISE_ALGN_PATH, 'rb') as _fh:
    accession_arr = np.load(_fh)
    cdd_master_seq_algn_mat = np.load(_fh)

cdd_master_seq_algn_feat_df = pd.DataFrame(
    cdd_master_seq_algn_mat,
    index=accession_arr,
)

cdd_superfamily_value_counts = \
    cdd_superfamily_df['superfamily_accession'].value_counts()
selected_cdd_superfamilies = \
    list(cdd_superfamily_value_counts[cdd_superfamily_value_counts > MIN_SUPERFAMILY_SIZE].index)
selected_cdd_superfamily_df = \
    cdd_superfamily_df[cdd_superfamily_df['superfamily_accession'].isin(selected_cdd_superfamilies)]


def visualize_conserved_domains(
        method: str,
        kwargs: dict,
):

    plt.figure(figsize=(24, 16))

    if method == 'tsne':
        pca = PCA(**kwargs)
        tsne = TSNE(n_components=2)
        cdd_master_seq_algn_coord = \
            tsne.fit_transform(pca.fit_transform(cdd_master_seq_algn_feat_df))
    elif method == 'umap':
        # if init option is not specified, it's 'spectral' by default
        # will throw 'graph is not fully connected' warning and hang
        umap = UMAP(**kwargs, init='random')
        cdd_master_seq_algn_coord = \
            umap.fit_transform(cdd_master_seq_algn_feat_df)
    else:
        raise ValueError(f'Unknown visualization method: {method}')

    if method == 'tsne':
        fig_name = \
            f'pca_n_components_{kwargs["n_components"]:03d}_' \
            f'tsne_n_components_2.png'
    else:
        fig_name = \
            f'umap_n_neighbors_{kwargs["n_neighbors"]:06d}_' \
            f'umap_min_dist_{kwargs["min_dist"]:.3f}.png'
    fig_dir_path = CDD_TSNE_IMAGE_DIR_PATH \
        if method == 'tsne' else CDD_UMAP_IMAGE_DIR_PATH
    fig_path = os.path.join(fig_dir_path, fig_name)

    cdd_master_seq_algn_coord_df = pd.DataFrame(
        cdd_master_seq_algn_coord,
        index=cdd_master_seq_algn_feat_df.index,
        columns=['x', 'y'],
    )

    cdd_master_seq_algn_coord_df.index.name = 'accession'
    cdd_coord_superfamily_df = pd.merge(
        cdd_master_seq_algn_coord_df,
        selected_cdd_superfamily_df[['superfamily_accession']],
        how='outer',
        on=['accession'],
    )
    cdd_coord_superfamily_df = cdd_coord_superfamily_df.fillna('others')

    num_superfamilies = \
        len(cdd_coord_superfamily_df['superfamily_accession'].unique())
    palette = ['#DDDDDD', ] + sns.color_palette('hls', num_superfamilies - 1)

    sns.scatterplot(
        data=cdd_coord_superfamily_df,
        x='x', y='y',
        hue='superfamily_accession',
        palette=palette,
        legend='brief',
    )
    plt.legend(loc='upper right', ncol=3)
    plt_title = \
        f'{"t-SNE" if method == "tsne" else "UMAP"} visualization of ' \
        f'conserved domains, with bigger (>= {MIN_SUPERFAMILY_SIZE} ' \
        f'members) superfamilies colored'
    plt.title(plt_title)
    plt.tight_layout()
    plt.savefig(fig_path)
    plt.cla()
    print(f'finished {fig_path}.')

# TODO: add histogram for conserved domains over all (selected) e coli genomes
# note that 'hit types' should be specified
# - specific hits
# - specific/superfamily (superfamily or specific without superfamily) hits
# - all (specific + non-specific) hits
# make a table of cols = [specific, specific/superfamily, all]


# script (executable) conserved domain visualization
if __name__ == '__main__':

    # plot the 2D space of conserved domains with t-SNE
    PCA_N_COMPONENTS_RANGE = list(range(2, 51))
    for _pca_n_components in PCA_N_COMPONENTS_RANGE:
        visualize_conserved_domains(
            'tsne',
            {
                'n_components': _pca_n_components,
            }
        )

    # plot the 2D space of conserved domains with UMAP
    UMAP_N_NEIGHBORS_RANGE = [
        (2 ** _i) for _i in
        range(1, int(np.ceil(np.log2(len(cdd_master_seq_algn_feat_df) / 4))))
    ]
    UMAP_MIN_DIST_RANGE = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    for _umap_n_neighbors in UMAP_N_NEIGHBORS_RANGE:
        for _umap_min_dist in UMAP_MIN_DIST_RANGE:
            visualize_conserved_domains(
                'umap', {
                    'n_neighbors': _umap_n_neighbors,
                    'min_dist': _umap_min_dist,
                }
            )
