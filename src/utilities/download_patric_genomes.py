"""
File Name:          download_patric_genomes.py
Project:            bioseq-learning

File Description:

"""
import os
import logging
from multiprocessing import Pool
from subprocess import Popen, PIPE
from typing import Optional, List, Tuple

from tqdm import tqdm

_LOGGER = logging.getLogger(__name__)
PATRIC_GENOMES_FTP_URL = 'ftp://ftp.patricbrc.org/genomes'


def download_patric_genome(
        genome_dir_path: str,
        extensions: List[str],
        genome_id: Optional[str],
):
    # get the genome ID if not given
    genome_dir_path: str = os.path.join(genome_dir_path, '')
    genome_id: str = genome_id if genome_id else \
        os.path.basename(genome_dir_path[:-1])

    # download the genome sequences from PATRIC
    _wget_cmd = Popen(
        [
            'wget',
            f'{PATRIC_GENOMES_FTP_URL}/{genome_id}/',
            '--no-clobber',
            '--recursive',
            '--no-parent',
            '--no-directories',
            '--accept', f'{",".join(extensions)}',
            '--quiet',
            '--directory-prefix', genome_dir_path,
        ],
        stdout=PIPE,
        stderr=PIPE,
        stdin=PIPE,
    )
    _wget_cmd.communicate()


def __download_patric_genome(arg):
    return download_patric_genome(*arg)


def download_patric_genomes(
        genome_parent_dir_path: str,
        genome_id_list: List[str],
        extensions: List[str],
        num_workers: int = 1,
):
    if not os.path.exists(genome_parent_dir_path):
        os.mkdir(genome_parent_dir_path)

    # clamp the number of workers between range [1, number of CPU cores]
    num_workers: int = max(1, min(num_workers, os.cpu_count()))

    # get all the arguments for genome downloading
    _genome_parent_dir_path = os.path.abspath(genome_parent_dir_path)
    _download_patric_genome_arguments: \
        List[Tuple[str, List[str], Optional[str]]] = []
    for _genome_id in genome_id_list:
        _genome_dir_path: str = \
            os.path.join(_genome_parent_dir_path, _genome_id)
        _download_patric_genome_arguments.append(
            (_genome_dir_path, extensions, _genome_id)
        )

    # parallel genome processing with pool
    with Pool(num_workers) as _pool:
        for _ in tqdm(
                _pool.imap_unordered(
                    __download_patric_genome,
                    _download_patric_genome_arguments,
                ),
                ncols=80,
                smoothing=0.1,
                total=len(_download_patric_genome_arguments)):
            pass