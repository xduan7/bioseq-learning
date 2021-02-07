#!/bin/bash

# change working directory to the directory of this script
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"

# process genomes in parallel
export PYTHONPATH=${PROJECT_DIR}:$PYTHONPATH
python ${PROJECT_DIR}/src/datasets/process_patric_faa_genomes.py \
  -i ${PROJECT_DIR}/data/raw/genomes/reference_or_representative_bacteria \
  -o ${PROJECT_DIR}/data/interim/genomes/reference_or_representative_bacteria \
  -w 160 \
  2> ${PROJECT_DIR}/logs/process_patric_ref_or_rep_bacteria_faa_genomes.txt
