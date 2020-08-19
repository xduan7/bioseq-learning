#!/bin/bash

# change working directory to the directory of this script
cd "${0%/*}"

# download NCBI conserved domain database and the internal tracking information
mkdir -p ../data/raw/CDD
wget -nc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz \
  -P ../data/raw/CDD/
wget -nc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt \
  -P ../data/raw/CDD/

# unpack and move the CDD files to interim
mkdir -p ../data/interim/CDD
tar -zxf ../data/raw/CDD/Cdd_LE.tar.gz \
  --directory ../data/interim/CDD
cp ../data/raw/CDD/cdtrack.txt ../data/interim/CDD/cd_track.txt
