#!/usr/bin/env bash


# Test CAVA suite. After having installed CAVA, execute this script as:
#
#   bash test.sh

set -x
set -e 
set -o pipefail

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
    samtools faidx data/tmp.GRCh38.fa
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
python3 bin/EnsemblDB.py -e 101 -o ENST101_test
python3 bin/EnsemblDB.py -e 101 -o ENST101_test_small -i test/CustomTX.txt
python3 bin/EnsemblDB.py -e 75 -o ENST75_test
python3 bin/EnsemblDB.py -e 75 -o ENST75_test_small -i test/CustomTX.txt


# Annotate
# build 38 test full
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp1 -c test/CAVA_config_1.txt
rm test/tmp1.log
# build 38 test small
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp2 -c test/CAVA_config_2.txt
rm test/tmp2.log
# build 37 test full
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp3 -c test/CAVA_config_3.txt
rm test/tmp3.log
# build 37 test small
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp4 -c test/CAVA_config_4.txt
rm test/tmp4.log

# ########################################################################################
# Test RefSeq
# ########################################################################################

# Prepare database
python3 bin/RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq
python3 bin/RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq_small -i test/CustomTX2.txt

# Annotate
# build 38 test full
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmpA -c test/CAVA_config_5.txt
rm test/tmpA.log
# build 38 test small
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmpB -c test/CAVA_config_6.txt
rm test/tmpB.log
# build 37 test full
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmpC -c test/CAVA_config_7.txt
rm test/tmpC.log
# build 37 test small
python3 bin/exit
CAVA.py -i test/test.GRCh37.vcf -o test/tmpD -c test/CAVA_config_8.txt
rm test/tmpD.log

# ########################################################################################
# Test MANE
# ########################################################################################

# Prepare database
python3 bin/MANE.py -e 0.91

# Annotate
# build 38 MANE ENST
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp_a -c test/CAVA_config_9.txt
rm test/tmp_a.log
# build 38 MANE REFSEQ
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp_b -c test/CAVA_config_10.txt
rm test/tmp_b.log
# build 37 MANE converted ENST
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp_c -c test/CAVA_config_11.txt
rm test/tmp_c.log
# build 37 MANE Converted hg19
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp_d -c test/CAVA_config_12.txt
rm test/tmp_d.log

echo "Running unit test for HGVSP"
python3 -m unittest test/test_csn.py

# Test: verify correct results
# =========
# Before comparing outputs, we have to remove the header from the VCFs, otherwise the md5sums won't match
for x in test/tmp*.vcf
do
  grep -v "#" $x > tmp.txt
  #cp tmp.txt ${x}.expected
  mv tmp.txt ${x}
done

md5sum --check test/hashes.txt

