#!/usr/bin/env bash


# Test CAVA suite. After having installed CAVA, execute this script as:
#
#   bash test.sh

set -x
set -e 
set -o pipefail

echo "Running unit test for HGVSP"
python3 -m unittest test/test_end2end.py

echo "Running unit tests for Variant"
python3 -m unittest test/test_csn.py

# Set up
#
# Download common variants to test 1% of all common variants as a robustness test.
# 
#if [ ! -f test/test.input.vcf.gz ]
#then
#    curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz \
#        | gunzip \
#        | awk '$1 ~ "^#" || NR % 100 == 0' \
#        | bgzip > test/test.input.vcf.gz
#fi
#
# Download reference genomes
#
if [ ! -f data/tmp.GRCh38.fa ]
then
    curl --cipher 'DEFAULT:!DH'  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > data/tmp.GRCh38.fa
    samtools faidx data/tmp.GRCh36.fa
fi

if [ ! -f data/tmp.hg19.fa ]
then
    curl --cipher 'DEFAULT:!DH' https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip > data/tmp.hg19.fa
    samtools faidx data/tmp.hg19.fa
fi

# ########################################################################################
# Test ENSEMBL
# ########################################################################################

# Prepare database
python3 EnsemblDB.py -e 101 -o ENST101_test
python3 EnsemblDB.py -e 101 -o ENST101_test_small -i test/CustomTX.txt
python3 EnsemblDB.py -e 75 -o ENST75_test
python3 EnsemblDB.py -e 75 -o ENST75_test_small -i test/CustomTX.txt

# ########################################################################################
# Test RefSeq
# ########################################################################################

# Prepare database
python3 RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq
python3 RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq_small -i test/CustomTX2.txt

# ########################################################################################
# Test MANE
# ########################################################################################

# Prepare database
python3 MANE.py -e 0.91
