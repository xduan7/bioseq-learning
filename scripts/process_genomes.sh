#!/bin/bash

# change working directory to the directory of this script
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"

# process genomes in parallel
export PYTHONPATH=${PROJECT_DIR}:$PYTHONPATH
nohup python ${PROJECT_DIR}/src/datasets/process_genomes.py \
  -i ${PROJECT_DIR}/data/raw/genomes/escherichia_coli \
  -o ${PROJECT_DIR}/data/interim/genomes/escherichia_coli \
  -t 240 \
  >> ${PROJECT_DIR}/data/interim/genomes/process_genomes.txt &
