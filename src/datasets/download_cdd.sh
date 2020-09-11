#!/bin/bash

# get the directory to project
PROJECT_DIR="$(dirname "$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")")"

# download NCBI conserved domain database and the internal tracking information
RAW_CDD_PATH=${PROJECT_DIR}/data/raw/CDD
mkdir -p ${RAW_CDD_PATH}
wget -nc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz \
  -P ${RAW_CDD_PATH}
wget -nc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt \
  -P ${RAW_CDD_PATH}

# unpack and move the CDD files to interim
INTERIM_CDD_PATH=${PROJECT_DIR}/data/interim/CDD
mkdir -p ${INTERIM_CDD_PATH}
tar -zxf ${RAW_CDD_PATH}/Cdd_LE.tar.gz \
  --directory ${INTERIM_CDD_PATH}
cp ${RAW_CDD_PATH}/cdtrack.txt ${INTERIM_CDD_PATH}/cd_track.txt
