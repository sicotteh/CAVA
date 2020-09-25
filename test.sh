#!/usr/bin/env bash

# Test CAVA suite. After having installed CAVA, execute this script as:
#
#   bash test.sh

set -e 
set -o pipefail

# Set up
if [ ! -f test/tmp.GRCh38.fa ]
then
    curl --cipher 'DEFAULT:!DH' http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > test/tmp.GRCh38.fa
    samtools faidx test/tmp.GRCh38.fa
fi

if [ ! -f test/tmp.hg19.fa ]
then
    curl --cipher 'DEFAULT:!DH' https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip > test/tmp.hg19.fa
    samtools faidx test/tmp.hg19.fa
fi

# Install program
python3 setup.py install

# Test: prepare database
python3 bin/EnsemblDB.py -e 101 -o ENST101_test
python3 bin/EnsemblDB.py -e 101 -o ENST101_test_small -i test/CustomTX.txt
python3 bin/EnsemblDB.py -e 75 -o ENST75_test
python3 bin/EnsemblDB.py -e 75 -o ENST75_test_small -i test/CustomTX.txt


# Test: annotate
# build 38 test full
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp1 -c test/CAVA_config_1.txt
# build 38 test small
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp2 -c test/CAVA_config_2.txt
# build 37 test full
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp3 -c test/CAVA_config_3.txt
# build 37 test small
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp4 -c test/CAVA_config_4.txt

tail -n2 test/tmp*vcf # None should have '.' in CSN

# Tear down
# =========

#rm test/tmp.*
