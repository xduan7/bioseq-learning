"""
File Name:          create_genome_example.py
Project:            bioseq-learning

File Description:

    This file is a script to create genome examples, which are sliced
    contigs of PATRIC genome

    example:
    $python create_genome_example.py \
        -i ./562.2283 \
        -s 128 -e 131 \
        -o ../../tests/examples/genomes

"""
import os
import argparse

import pandas as pd
from Bio import SeqIO


parser = argparse.ArgumentParser(
    description='create a small example of PATRIC genomes for testing')
parser.add_argument(
    '-i', '--genome_dir_path', type=str,
    help='directory path to a PATRIC genome',
)
parser.add_argument(
    '-s', '--start', type=int,
    help='first index for contigs (1-based numbering)',
)
parser.add_argument(
    '-e', '--end', type=int,
    help='last index for contigs (1-based numbering, inclusive)',
)
parser.add_argument(
    '-o', '--output_dir_path', type=str,
    help='directory path to a PATRIC genome',
)
args = parser.parse_args()

genome_dir_path: str = os.path.join(args.genome_dir_path, '')
print(genome_dir_path)
genome_id: str = os.path.basename(genome_dir_path[:-1])
print(genome_id)

contig_seq_path: str = f'./{args.genome_dir_path}/{genome_id}.fna'
feature_path: str = \
    f'./{args.genome_dir_path}/{genome_id}.PATRIC.features.tab'

_start, _end = args.start - 1, args.end

_output_dir = os.path.join(args.output_dir_path, f'{genome_id}')
if not os.path.isdir(_output_dir):
    os.makedirs(_output_dir)

# load, slice and store the contigs
contig_seq_records = list(SeqIO.parse(contig_seq_path, 'fasta'))
_example_contig_seq_records = \
    [contig_seq_records[i].id for i in range(_start, _end)]
with open(os.path.join(_output_dir, f'{genome_id}.fna'), 'w+') as _fh:
    SeqIO.write(contig_seq_records[_start: _end], _fh, 'fasta')

# load, slice and store the features
feature_df = pd.read_table(feature_path)
_example_feature_df = feature_df[feature_df['accession'].isin([
    contig_seq_records[_i].id for _i in range(_start, _end)
])]
_example_feature_df.to_csv(
    os.path.join(_output_dir, f'{genome_id}.PATRIC.features.tab'),
    index=None,
    sep='\t',
)
