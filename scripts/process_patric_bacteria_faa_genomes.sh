#!/bin/bash

# change working directory to the directory of this script
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"

# process genomes in parallel
export PYTHONPATH=${PROJECT_DIR}:$PYTHONPATH
python ${PROJECT_DIR}/src/datasets/process_patric_faa_genomes.py \
  -i ${PROJECT_DIR}/data/raw/genomes/bacteria \
  -o ${PROJECT_DIR}/data/interim/genomes/bacteria \
  -w 160 \
  2> ${PROJECT_DIR}/logs/process_patric_faa_genomes.txt &
