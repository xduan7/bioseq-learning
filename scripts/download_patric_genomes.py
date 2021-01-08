"""
File Name:          download_patric_genomes.py.py
Project:            bioseq-learning

File Description:

    Download all the genomes from an automatically generated PATRIC table,
    with a list of accepted extensions (e.g. features.tab, PATRIC.faa, etc.)
    in case there are too many items and PATRIC refuses to initialize the
    direct download.

"""
import os
import sys
import argparse
import pandas as pd


if __name__ == '__main__':
    # usage example:
    # $python download_patric_genomes.py \
    #      -i ../data/raw/genomes/bacteria.csv \
    #      -o ../data/raw/genomes \
    #      -e features.tab faa fna \
    #      -w 80

    parser = argparse.ArgumentParser(
        description='download PATRIC genomes in parallel')

    parser.add_argument(
        '-i', '--input_genome_table', type=str, required=True,
        help='path to the table (*.csv) file with \'Genome ID\' column',
    )
    parser.add_argument(
        '-o', '--output_genome_parent_dir_path', type=str, required=True,
        help='parent directory path for genomes, the files for each genome '
             'will be stored in the subdirectory named after its genome ID',
    )
    parser.add_argument(
        '-e', '--extensions', nargs='*', type=str, default=['features.tab'],
        help='list of extensions to download, e.g. fna, PATRIC.faa, etc.',
    )
    parser.add_argument(
        '-w', '--num_workers', type=int, default=1,
        help='number of workers for parallelization',
    )
    args = parser.parse_args()

    _genome_df = pd.read_csv(
        args.input_genome_table,
        index_col=False,
        usecols=['Genome ID'],
        dtype={'Genome ID': str},
    )
    _genome_id_list = _genome_df['Genome ID'].tolist()

    # append the project source path to PYTHONPATH
    _project_src_dir = os.path.abspath(os.path.join(__file__, '../../'))
    sys.path.append(_project_src_dir)
    from src.datasets.download_patric_genomes import download_patric_genomes

    download_patric_genomes(
        genome_parent_dir_path=args.output_genome_parent_dir_path,
        genome_id_list=_genome_id_list,
        extensions=args.extensions,
        num_workers=args.num_workers,
    )
