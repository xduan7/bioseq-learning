#!/bin/bash

# get the paths to blast binary, the project, and the data directories
RPSTBLASTN="rpstblastn"
if ! command -v ${RPSTBLASTN} &> /dev/null; then
  printf "\e[1;31m[ERROR]\e[0m Could not find NCBI-Blast+ binary ${RPSTBLASTN}. Make sure that Blast binaries are included in the PATHS environment variable.\n"
  exit 1
fi
BLAST_BIN_DIR="$(dirname "$(which ${RPSTBLASTN})")"
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"
RAW_DATA_DIR=${PROJECT_DIR}/data/raw
INTERIM_DATA_DIR=${PROJECT_DIR}/data/interim
PROCESSED_DATA_DIR=${PROJECT_DIR}/data/processed

# download rpsbproc from NCBI FTP server
# note that this version is only executable on linux machines, and the source code is not compilable on Mac
RAW_RPSBPROC_DIR=${RAW_DATA_DIR}/RpsbProc
INTERIM_RPSBPROC_DIR=${INTERIM_DATA_DIR}/RpsbProc
mkdir -p "${RAW_RPSBPROC_DIR}"
mkdir -p "${INTERIM_RPSBPROC_DIR}"
wget -nc -q ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz -P "${RAW_RPSBPROC_DIR}"
tar -zxf "${RAW_RPSBPROC_DIR}"/RpsbProc-x64-linux.tar.gz --directory "${INTERIM_RPSBPROC_DIR}"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cp "${INTERIM_RPSBPROC_DIR}"/RpsbProc-x64-linux/rpsbproc "${BLAST_BIN_DIR}"
else
  printf "\e[1;31m[WARNING]\e[0m Downloaded RpsbProc from ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/RpsbProc-x64-linux.tar.gz is only executable on linux systems.\n"
fi

# download and unpack CDD (conserved domain database)
RAW_CDD_DIR=${RAW_DATA_DIR}/CDD
INTERIM_CDD_DIR=${INTERIM_DATA_DIR}/CDD
mkdir -p "${RAW_CDD_DIR}"
mkdir -p "${INTERIM_CDD_DIR}"
wget -nc -q ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz -P "${RAW_CDD_DIR}"
tar -zxf "${RAW_CDD_DIR}"/Cdd_LE.tar.gz --directory "${INTERIM_CDD_DIR}"

# download metadata files for CDD (required for rpsbproc)
RAW_CDD_METADATA_DIR=${RAW_DATA_DIR}/CDD_metadata
INTERIM_CDD_METADATA_DIR=${INTERIM_DATA_DIR}/CDD_metadata
mkdir -p "${RAW_CDD_METADATA_DIR}"
mkdir -p "${INTERIM_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt -P "${RAW_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt -P "${RAW_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links -P "${RAW_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P "${RAW_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz -P "${RAW_CDD_METADATA_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz -P "${RAW_CDD_METADATA_DIR}"
cp "${RAW_CDD_METADATA_DIR}"/* "${INTERIM_CDD_METADATA_DIR}"
gzip -d -k -f "${INTERIM_CDD_METADATA_DIR}"/*.gz
