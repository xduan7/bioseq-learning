"""
File Name:          conserved_domain_visualization.py
Project:            bioseq-learning

File Description:

"""
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from src.utilities import create_directory
from src import INTERIM_DATA_DIR_PATH, PROCESSED_DATA_DIR_PATH, DOC_DIR_PATH


MIN_SUPERFAMILY_SIZE = 128

PCA_N_COMPONENTS_RANGE = list(range(2, 51))
UMAP_N_NEIGHBORS_RANGE = [
    (2 ** _i) for _i in range(1, int(np.ceil(np.log2(55000 / 4))))]
UMAP_MIN_DIST_RANGE = list(np.arange(0.09, 1, 0.1))


CDD_SUPERFAMILY_PATH = os.path.join(
    INTERIM_DATA_DIR_PATH, 'CDD_metadata/family_superfamily_links')
CDD_MASTER_PROCESSED_ALGN_FEAT_PATH = os.path.join(
    PROCESSED_DATA_DIR_PATH, 'CDD_alignment/cdd_master_alignment_feature.csv')

CDD_TSNE_IMAGE_DIR_PATH = os.path.join(DOC_DIR_PATH, 'images/cdd_tsne')
CDD_UMAP_IMAGE_DIR_PATH = os.path.join(DOC_DIR_PATH, 'images/cdd_umap')


create_directory(CDD_TSNE_IMAGE_DIR_PATH)
create_directory(CDD_UMAP_IMAGE_DIR_PATH)

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

cdd_master_seq_algn_feat_df = pd.read_csv(
    CDD_MASTER_PROCESSED_ALGN_FEAT_PATH,
    index_col=0,
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
        cdd_master_seq_algn_coord = tsne.fit_transform(pca.fit_transform(
            cdd_master_seq_algn_feat_df))
    else:
        umap = UMAP(**kwargs)
        cdd_master_seq_algn_coord = umap.fit_transform(cdd_master_seq_algn_feat_df)

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
    plt.title(f'{method} ({kwargs}) visualization of conserved domains')
    plt.tight_layout()
    plt.savefig(os.path.join(
        CDD_TSNE_IMAGE_DIR_PATH if method == 'tsne'
        else CDD_UMAP_IMAGE_DIR_PATH,
        f'{str(kwargs).replace(" ", "_").replace(":", "_")}.png',
    ))
    plt.cla()


for _pca_n_components in PCA_N_COMPONENTS_RANGE:
    visualize_conserved_domains('pca', {'n_components': _pca_n_components})

for _umap_n_neighbors in UMAP_N_NEIGHBORS_RANGE:
    for _umap_min_dist in UMAP_MIN_DIST_RANGE:
        visualize_conserved_domains(
            'tsne', {
                'n_neighbors': _umap_n_neighbors,
                'min_dist': _umap_min_dist,
            }
        )

