"""
File Name:          visualize_reference_or_representative_bacteria.py
Project:            bioseq-learning

File Description:

Visualization of Reference or Representative Bacteria
- [ ] organism distribution (pie chart)
- [ ] organism distribution (pie chart)
- [ ] number of contigs per genome distribution (histogram)
- [ ] number of conserved domains per coding region
- [ ] conserved domain frequency (maybe a pie chart)

"""
import os
import pickle

import taxonomy
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List, Optional, Tuple

from src import NCBI_TAX_DIR_PATH
from src.datasets.genome_domain_dataset import (
    Annotation,
    ContigWithConservedDomains,
    Organism,
)


ANNOTATION = Annotation.PATRIC
COLORS = plt.cm.Pastel1.colors

REF_OR_REP_BACTERIA_CONTIGS_WITH_CDS_FILE_PATH = \
    f'../data/processed/genomes/reference_or_representative_' \
    f'bacteria_contigs_with_conserved_domains.{ANNOTATION.value}.pickle'

with open(REF_OR_REP_BACTERIA_CONTIGS_WITH_CDS_FILE_PATH, 'rb') as __fh:
    contigs_with_cds: List[ContigWithConservedDomains] = pickle.load(__fh)

# @dataclass
# class ContigWithConservedDomains:
#     """Data class for a genome contig, annotated with features (PATRIC or
#     RefSeq) and the corresponding conserved domains.
#     """
#     genome_id: str
#     genome_name: Optional[str]
#     organism: Optional[Organism]
#     annotation: Annotation
#     contig_accession: str
#     contig_feature_df: pd.DataFrame
#     contig_conserved_domain_df: pd.DataFrame
#     contig_feature_csv_file_path: str
#     contig_conserved_domain_csv_file_path: str

#

_contig_info_list: \
    List[Tuple[str, Optional[str], Optional[Organism], str, str]] = []
for __contig_with_cds in contigs_with_cds:
    __contig_info = (
        __contig_with_cds.genome_id,
        __contig_with_cds.genome_name,
        __contig_with_cds.organism.value,
        __contig_with_cds.ncbi_taxon_id,
        __contig_with_cds.contig_accession
    )
    _contig_info_list.append(__contig_info)
contig_info_df = pd.DataFrame(
    _contig_info_list,
    columns=[
        'genome_id',
        'genome_name',
        'organism',
        'ncbi_taxon_id',
        'contig_accession',
    ],
)

__print_msg = \
    f'There are {len(contig_info_df["genome_id"].unique())} genomes with ' \
    f'{len(contig_info_df)} contigs from all the reference or ' \
    f'representative bacteria with {ANNOTATION.value} annotation.'
print(__print_msg)


# visualization of organism distribution pie chart
_fig1_path = \
    f'../docs/images/reference_or_representative_bacteria/' \
    f'organism_pie.{ANNOTATION.value}.png'
if not os.path.exists(_fig1_path):
    print(f'generating {ANNOTATION.value} organism distribution pie chart ...')
    _fig1, (_fig1_ax1, _fig1_ax2) = plt.subplots(1, 2, figsize=(48, 24))
    _fig1.tight_layout(pad=8)
    _fig1_ax1.pie(
        contig_info_df.drop_duplicates(subset=['genome_id'])['organism'].value_counts().values,
        labels=contig_info_df.drop_duplicates(subset=['genome_id'])['organism'].value_counts().index,
        colors=COLORS,
        autopct='%1.1f%%',
        pctdistance=0.75,
        labeldistance=0.8,
    )
    _fig1_ax1.axis('equal')
    _fig1_ax1.set_title(
        f'Organism Percentage by '
        f'{len(contig_info_df["genome_id"].unique())} Genomes'
    )
    _fig1_ax2.pie(
        contig_info_df['organism'].value_counts().values,
        labels=contig_info_df['organism'].value_counts().index,
        colors=COLORS,
        autopct='%1.1f%%',
        pctdistance=0.75,
        labeldistance=0.8,
    )
    _fig1_ax2.axis('equal')
    _fig1_ax2.set_title(f'Organism Percentage by {len(contig_info_df)} Contigs')

    os.makedirs(os.path.dirname(_fig1_path), exist_ok=True)
    plt.savefig(_fig1_path)


# visualization of organism distribution sun burst
_fig2_path = \
    f'../docs/images/reference_or_representative_bacteria/' \
    f'organism_sunburst.{ANNOTATION.value}.html'
if not os.path.exists(_fig2_path):
    print(f'generating {ANNOTATION.value} organism distribution sun burst ...')
    tax = taxonomy.Taxonomy.from_ncbi(
        nodes_path=os.path.join(NCBI_TAX_DIR_PATH, 'nodes.dmp'),
        names_path=os.path.join(NCBI_TAX_DIR_PATH, 'names.dmp'),
    )
    taxon_ids = contig_info_df['ncbi_taxon_id'].drop_duplicates().tolist()
    taxon_info = []
    taxon_ranks = [
        'species',
        'genus',
        # 'tribe',
        # 'subfamily',
        'family',
        'order',
        # 'clade',
        # 'subphylum',
        'phylum',
        # 'kingdom',
    ]

    for __ncbi_taxon_id in tqdm(taxon_ids):
        __ncbi_taxon_id = str(__ncbi_taxon_id)
        __taxon_info = [__ncbi_taxon_id, ]
        for __rank in taxon_ranks:
            __parent = tax.parent(__ncbi_taxon_id, __rank)
            __taxon_info.append(__parent.name if __parent else '')
        taxon_info.append(__taxon_info)

    tax_info_df = pd.DataFrame(
        taxon_info,
        columns=['ncbi_taxon_id'] + taxon_ranks,
    )
    fig2 = px.sunburst(
        tax_info_df,
        path=reversed(taxon_ranks),
    )
    fig2.write_html(_fig2_path)