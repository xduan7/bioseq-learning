"""
File Name:          download_patric_genomes.py
Project:            bioseq-learning

File Description:

    This file implements the download methods for genomes on PTRAIC FTP site.

"""
import os
import logging
import subprocess
from multiprocessing import Pool
from typing import Optional, List, Tuple

from tqdm import tqdm

_LOGGER = logging.getLogger(__name__)
MAX_NUM_ATTEMPTS = 10
RECOMMENDED_NUM_WORKERS = 5
PATRIC_GENOMES_FTP_URL = 'ftp://ftp.patricbrc.org/genomes'


def download_patric_genome(
        genome_dir_path: str,
        extensions: List[str],
        genome_id: Optional[str],
):
    """download a single genome from PATRIC FTP site

    :param genome_dir_path: path to the directory that stores the genome
    sequence files (e.g. faa, fna, etc.)
    :type genome_dir_path: str
    :param extensions: file extensions that will be downloaded (e.g.
    PATRIC.faa, features.tab, fna, etc.)
    :type extensions: List[str]
    :param genome_id: PATRIC genome ID; if not given, it will be infereed
    from the base path name of the genome_dir_path argument
    :type genome_id: Optional[str]
    :return: None
    :rtype: None
    """
    # get the genome ID if not given
    genome_dir_path: str = os.path.join(genome_dir_path, '')
    genome_id: str = genome_id if genome_id else \
        os.path.basename(genome_dir_path[:-1])

    # download the genome sequences from PATRIC
    _wget_cmd_list = [
        'wget',
        '--no-clobber',
        '--recursive',
        '--no-parent',
        '--no-directories',
        '--accept', f'{",".join(extensions)}',
        # '--quiet',  # not set to get the error message
        '--directory-prefix', genome_dir_path,
        f'{PATRIC_GENOMES_FTP_URL}/{genome_id}/',
    ]
    # _wget_cmd_str = ' '.join(_wget_cmd_list)
    _num_attempt = 0
    while _num_attempt < MAX_NUM_ATTEMPTS:
        _wget_cmd = subprocess.Popen(
            _wget_cmd_list,
            text=True,
            # shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        _wget_cmd.wait()
        _wget_stdout, _wget_stderr = _wget_cmd.communicate()
        if _wget_cmd.returncode != 0:
            _LOGGER.warning(f'wget error: {_wget_stderr}')
            _num_attempt += 1
        else:
            return

    _warning_msg = \
        f'Downloading failed for {genome_id}: ' \
        f'maximum number of attempts exceeded.'
    _LOGGER.warning(_warning_msg)


def __download_patric_genome(arg):
    return download_patric_genome(*arg)


def download_patric_genomes(
        genome_parent_dir_path: str,
        genome_id_list: List[str],
        extensions: List[str],
        num_workers: int = 1,
):
    """download multiple patric genomes in parallel

    :param genome_parent_dir_path: parent directory of all the downloaded
    genomes; each genome will be stored inside the child directory named
    after their PATRIC genome ID
    :type genome_parent_dir_path: str
    :param genome_id_list: list of PATRIC genome IDS to be downloaded
    :type genome_id_list: List[str]
    :param extensions: file extensions that will be downloaded (e.g.
    PATRIC.faa, features.tab, fna, etc.)
    :type extensions: List[str]
    :param num_workers: number of workers for parallel downloading; note
    that too many workers might result in login error from the FTP server
    :type num_workers: int
    :return: None
    :rtype: None
    """
    if not os.path.exists(genome_parent_dir_path):
        os.mkdir(genome_parent_dir_path)

    # clamp the number of workers between range [1, number of CPU cores]
    num_workers: int = max(1, min(num_workers, os.cpu_count()))
    if num_workers > RECOMMENDED_NUM_WORKERS:
        _warning_msg = \
            'Too many downloading workers, might raise the ' \
            '\'incorrect login\' error from PATRIC FTP server.'
        _LOGGER.warning(_warning_msg)

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
