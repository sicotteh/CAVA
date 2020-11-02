#!/usr/bin/env bash

# Test CAVA suite. After having installed CAVA, execute this script as:
#
#   bash test.sh

set -x
set -e 
set -o pipefail

# Set up
# changed to test 1% of all common variants.
# 
if [ ! -f test/test.input.vcf.gz ]
then
    curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz \
        | gunzip \
        | awk '$1 ~ "^#" || NR % 100 == 0' \
        | bgzip > test/test.input.vcf.gz
fi

#if [ ! -f test/tmp.hg38.fa ]
#then
#    curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > test/tmp.hg38.fa
#    samtools faidx test/tmp.hg38.fa
#fi

# Test: prepare database
#./ensembl_db -e 75 -o test/tmp.db

# Lengths in bp not kb:
#gunzip -c test/tmp.db.gz | grep -F '+/919bp/1/918bp/305' > /dev/null

# Test: annotate
#CONFIG=test/CAVA_config.txt
CONFIG=/research/bsi/projects/staff_analysis/m037385/CAVA/cavaconfig.GRCh38.txt
echo "Testing Integration test list"
./cava -i test/integration_test_278.sorted.vcf -o test/test.out -c ${CONFIG} -t 8 -o integration_test_278.sorted.out
echo "Testing a lot of common variants"
./cava -i test/test.input.vcf.gz -o test/test.out -c ${CONFIG} -t 8



# Tear down
# =========

rm test/tmp.*
