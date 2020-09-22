#!/bin/bash

# change working directory to the directory of this script
PROJECT_DIR="$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")"

# align conserved domains (using master sequences)
export PYTHONPATH=${PROJECT_DIR}:$PYTHONPATH
nohup python ${PROJECT_DIR}/src/datasets/conserved_domain_alignment.py \
  >> ${PROJECT_DIR}/logs/align_conserved_domains.txt &
