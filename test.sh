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
if [ ! -f test/test.input.vcf.gz ]
then
    curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz \
        | gunzip \
        | awk '$1 ~ "^#" || NR % 100 == 0' \
        | bgzip > test/test.input.vcf.gz
fi
#
# Download reference genomes 
#
if [ ! -f test/tmp.GRCh38.fa ]
then
    curl  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > test/tmp.GRCh38.fa
    samtools faidx test/tmp.GRCh38.fa
fi

if [ ! -f test/tmp.hg19.fa ]
then
    curl  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gunzip > test/tmp.hg19.fa
    samtools faidx test/tmp.hg19.fa
fi


# Install program


ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
echo "ABSOLUTE_PATH=${ABSOLUTE_PATH}"
unset PYTHONPATH
export PYTHONPATH=${ABSOLUTE_PATH}/lib/python3.7/site-packages/
python3 setup.py build -f 
python3 setup.py install --prefix "${ABSOLUTE_PATH}"


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


# ########################################################################################
# Test RefSeq
# ########################################################################################

# Test: prepare database
python3 bin/RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq
python3 bin/RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq_small -i test/CustomTX2.txt
# The script to prepare a mapping between NM and NP is really slow (1 per second), so we're only
# running it for the small database (CAVA_config_7.txt)
bash NM2NP.bash data/RefSeq_small.gz data/RefSeq_small.nmnp.txt test/TMP

# Test: annotate
# build 38 test full
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmpA -c test/CAVA_config_5.txt
# build 38 test small
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmpB -c test/CAVA_config_6.txt
# build 37 test full
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmpC -c test/CAVA_config_7.txt
# build 37 test small
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmpD -c test/CAVA_config_8.txt


# ########################################################################################
# Test MANE
# ########################################################################################

python3 bin/MANE.py -e 0.91
# build 38 MANE ENST
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp_a -c test/CAVA_config_9.txt
# build 38 MANE REFSEQ
python3 bin/CAVA.py -i test/test.GRCh38.vcf -o test/tmp_b -c test/CAVA_config_10.txt
# build 37 MANE converted ENST
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp_c -c test/CAVA_config_11.txt
# build 37 MANE Converted hg19
python3 bin/CAVA.py -i test/test.GRCh37.vcf -o test/tmp_d -c test/CAVA_config_12.txt


echo "Running unit test for HGVSP"
PWD=`pwd`
cd cava_/
python3 -m unittest test_csn.py
cd ${PWD}

echo "Testing edge cases"
# Test: 
echo "Testing edge cases"

CONFIG=cavaconfig.GRCh38.txt
./cava -i test/integration_test_278.sorted.vcf -o test/integration_test_278.sorted.out -c ${CONFIG} -t 8


# Test: verify correct results
# =========
# Before comparing outputs, we have to remove the date from the VCFs, otherwise the md5sums won't match
for x in test/tmp*.vcf
do
  grep -v fileDate $x > tmp.txt
  #cp tmp.txt ${x}.expected
  mv tmp.txt ${x}
done


md5sum --check test/hashes.txt



echo "Testing a lot of common variants for runtime stability"
./cava -i test/test.input.vcf.gz -o test/test.common.out -c ${CONFIG} -t 8

