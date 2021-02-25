#!/bin/bash

# get the paths to the project and the data directories
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"
RAW_DATA_DIR=${PROJECT_DIR}/data/raw
INTERIM_DATA_DIR=${PROJECT_DIR}/data/interim

RAW_NCBI_TAX_DIR=${RAW_DATA_DIR}/NCBI_taxonomy
INTERIM_NCBI_TAX_DIR=${INTERIM_DATA_DIR}/NCBI_taxonomy
mkdir -p "${RAW_NCBI_TAX_DIR}"
mkdir -p "${INTERIM_NCBI_TAX_DIR}"
wget -nc -q https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -P "${RAW_NCBI_TAX_DIR}"
tar -zxf "${RAW_NCBI_TAX_DIR}"/new_taxdump.tar.gz --directory "${INTERIM_NCBI_TAX_DIR}"
