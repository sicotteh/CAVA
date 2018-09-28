#!/usr/bin/env bash

# Test CAVA suite. After having installed CAVA, execute this script as:
#
#   bash test.sh

set -e 
set -o pipefail

# Set up

if [ ! -f test/tmp.input.vcf.gz ]
then
    curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz \
        | gunzip \
        | awk '$1 ~ "^#" || NR % 1000 == 0' \
        | gzip > test/tmp.input.vcf.gz
fi

if [ ! -f test/tmp.hg38.fa ]
then
    curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > test/tmp.hg38.fa
    samtools faidx test/tmp.hg38.fa
fi

# Test: prepare database
./ensembl_db -e 75 -o test/tmp.db

# Lengths in bp not kb:
gunzip -c test/tmp.db.gz | grep -F '+/919bp/1/918bp/305' > /dev/null

# Test: annotate
./cava -i test/tmp.input.vcf.gz -o test/tmp. -c test/CAVA_config.txt -t 8

# Tear down
# =========

rm test/tmp.*
