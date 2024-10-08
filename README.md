 

* [CAVA README](#cava-readme)
    * [1 INTRODUCTION](#1-introduction)
    * [2 PUBLICATION](#2-publication)
    * [3 DEPENDENCIES](#3-dependencies)
    * [4 INSTALLATION ON LINUX OR MAC](#4-installation-on-linux-or-mac)
    * [5 RUNNING CAVA](#5-running-cava)
    * [6 LICENCE](#6-licence)

CAVA README
==================

1 INTRODUCTION
--------------

CAVA (Clinical Annotation of VAriants) is a lightweight, fast, flexible 
and easy-to-use Next Generation Sequencing (NGS) variant annotation tool. 
It implements a clinical sequencing nomenclature (CSN), a fixed variant 
annotation consistent with the principles of the Human Genome Variation 
Society (HGVS) guidelines, optimised for automated clinical variant 
annotation of NGS data. 

Since 2017, CAVA has been maintained by a group of bioinformaticians involved 
in both research and clinical genomics, adding several key functionalities, enhancements, 
and general support along the way.

2 PUBLICATION
-------------

If you use CAVA, please cite:

Márton Münz, Elise Ruark, Anthony Renwick, Emma Ramsay, Matthew Clarke, 
Shazia Mahamdallie, Victoria Cloke, Sheila Seal, Ann Strydom, 
Gerton Lunter, Nazneen Rahman. CSN and CAVA: variant annotation tools 
for rapid, robust next-generation sequencing analysis in the clinical 
setting. Genome Medicine 7:76, doi:10.1186/s13073-015-0195-6 (2015).

Maybe some day, we'll get around to publishing what we've done to rescue this abandoned project.

3 DEPENDENCIES
--------------

To install and run CAVA you will need the following dependencies installed:
- Python 3
- GCC and GNU make
- virtualenv

It just makes sense to keep things in a virtualenv. Here's how you do it if you are
unfamiliar.

```bash 
pip install virtualenv
virtualenv cava
source cava/bin/activate
```
At this point, your terminal should change to let you know you are in a virtual environment.

4 INSTALLATION ON LINUX OR MAC
------------------------------

```bash 
pip install cava

# - or -
git clone git@github.com:Steven-N-Hart/CAVA.git
# optional to checkout release
# e.g. git checkout v.1.2.4
cd CAVA
python setup.py install
```

If you get an error with pycurl, run the following command before running the setup.py
```bash
pip uninstall pycurl
export PYCURL_SSL_LIBRARY=nss
pip install --compile --install-option="--with-nss" --no-cache-dir pycurl
#OR (depending on the python version)
pip install --compile --global-option="--with-nss" --no-cache-dir pycurl  
```


5 RUNNING CAVA
--------------

Before using CAVA, you will need to create a config file. You have to provide two main components.
1) A fasta reference file and matching index file .fai
   If you do not have a file in your space.
   cd CAVA/cava/data
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 
   gunzip -c hg38.fa.gz > tmp.GRCh38.fa
   samtools faidx  tmp.GRCh38.fa

2) Create a database of transcripts for which to base your annotations from.
Details can be found in [this README](cava/ensembldb/README.md). In short, we recomend using MANE transcripts, 
so to get started, you would simply:
```bash
# Download GTF files for either RefSeq or ENSEMB
# Option 1: Use our script (after you install CAVA)
python3 MANE.py --no_hg19 -e 1.1 --outdir data


# Option 2: Download Manually (adjust version numbers to latest)
wget -O data/ENST.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.1/MANE.GRCh38.v1.1.ensembl_genomic.gtf.gz and
wget -O data/RefSeq.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.1/MANE.GRCh38.v1.1.refseq_genomic.gtf.gz

# Separate into ENST and NM Transcripts
zcat data/ENST.gtf.gz |cut -f9|cut -f4 -d' '|grep ENST|sed 's/;//;s/\"//g'|sort -u > data/ENST.txt
zcat data/RefSeq.gtf.gz |cut -f9|cut -f4 -d' '|grep "NM_"|sed 's/;//;s/\"//g'|sort -u > data/RefSeq.txt
```
Finally, create a config.txt files using the provided config_template.txt, but provide the location of the fasta reference and the ensembl transcript database you dowloaded. CAVA then can be run with the following simple command.

```bash
python3 CAVA.py -c config.txt -i input.vcf -o output
```

It requires three command line arguments: 
the name of the configuration file (-c), the name of the input file (-i) 
and the prefix of the output file name (-o). 

6 LICENCE
---------

CAVA is released under MIT licence (see the LICENCE file).

7 CHANGES HISTORY
---------
This version of CAVA includes the following changes (aside from bug fixes, especially for edge cases where multiple interpretations could apply)
- support refseq transcripts in addition to ensembl
- include new tags: CAVA_HGVSG, CAVA_HGVSC, CAVA_HGVSP to represent the current full HGVS nomenclature (G=Genomics, C=CDNA,P=Protein) for the HGVSC and HGVSP. We do not support the genomic tandem repeats for HGVSG nor imperfect repeats. These fields must be URL-decode (uudecode) because they include ';' encoded as %3B (';' is not a legal character in the VCF INFO field).
   -- This requires the addition of a files to support the protein information. Catalogs prior to version 2.03 will not be compatible because of that additional file.
- Include support for selenocysteine genes and alternate stop codon. Usually this is a conservative interpretation, often resulting in p.? when a novel stop codon might be discovered.
- Support for MANE transcripts (including transcripts with UTR's of length 0  which are part of MANE 1.0)
- Includes a lot of unit test for edge cases.
- Known Limitations: 
      -- large variants with coordinates outside transcript(we have partial support for those). 
      -- Novel Start Gain in 5'UTR are not reported (the impact of those putative changes is hard to predict.
      -- The %3B splitting may be hard when a variant has repeats and multiple transcripts.
- Optimized for speed via recoding and caching. From 20-500 times faster
- Support Selenocysteine genes in a conservative fashion. Even though we could get annotation for the location of the SECIS element from NCBI, we don't know the range of the effectiveness of the SECIS element. We know that the SECIS element no longer works at or past the original UAG Stop Codon.. so any new stop codon will be recoded as a U (selenocysteine) unless the protein becomes longer than the original one, so frameshifts or mutation/deletion of the original Stop codon will eventually find another stop codon that will not be recoded as U.
2.0.12 Changes
-Bug fixes for insertions/deletions that can be normalized right at the intron/exon junction edge. Were missing SO values and were sometimes called INT.
Introduced Initiator Gain (IG) (aka Start-Gain) features in 5'UTR for novel Start codons created upstream and in phase of the cannonical AUG. Annotated with CLASS=IG. Must update the impactdef tag in the config file to include an IMPACT for IG (currently IMPACT=3). Note that there is no SO equivalent for the IG tag.
- Better support for large deletions spanning intron/exon junction.



