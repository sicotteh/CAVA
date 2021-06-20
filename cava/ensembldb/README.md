# Preparing a transcript database

## Step 1. Identify a transcript set (optional, but recommended)

The best place to get a list of transcript to use is from the [MANE](ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/)
project. These are curated & versioned lists of what should be considered the "default" transcript for most genes.

In this case, we will use release 0.91, but other options are listed further down.

```
# Download GTF files for either RefSeq or ENSEMBLE
wget -O data/ENST.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.91/MANE.GRCh38.v0.91.select_ensembl_genomic.gtf.gz and
wget -O data/RefSeq.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.91/MANE.GRCh38.v0.91.select_refseq_genomic.gtf.gz

# Separate into ENST and NM Transcripts
zcat data/ENST.gtf.gz |cut -f9|cut -f4 -d' '|grep ENST|sed 's/;//;s/\"//g'|sort -u > data/ENST.txt
zcat data/RefSeq.gtf.gz |cut -f9|cut -f4 -d' '|grep "NM_"|sed 's/;//;s/\"//g'|sort -u > data/RefSeq.txt

# Look at the options and configure appropriately
python3 MANE.py -h

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -D OUTPUT_DIR, --outdir=OUTPUT_DIR
                        Output directory
  -e ENSEMBL, --ensembl=ENSEMBL
                        release version
  --no_hg19             Set this to skip hg19 builds

Example usage: MANE.py -e 0.91 -o mane_0.91
Note: by default, hg19 will be created using crossmap
Version: 1.3.3

```

## Option 2a. Make an Ensembl Database

To make an Ensemble database, use the following program:

```
Usage:

EnsemblDB.py

CAVA ensembl_db v1.3.0

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input filename (list of ENST IDs)
  -o OUTPUT, --output=OUTPUT
                        Output filename prefix
  -D OUTPUT_DIR, --outdir=OUTPUT_DIR
                        Output directory
  -e ENSEMBL, --ensembl=ENSEMBL
                        Ensembl release version.py
  --no_hg19             Set this to skip hg19 builds

Example usage: CAVA/ensembl_db -e 75 -o ensembl_db_75
Note: by default, hg19 will be created using crossmap
Version: 1.3.0
```

Note that by default, if you specify a version > 75, the coordinates will be in GRCh38. Additionally, this script will
automatically create an hg19 coordinate version unless you tell it not to. Using our example from above, the command

```
python3 EnsemblDB.py -e 101 -o ENST75 -D data -i data/ENST.txt
``` 

will create two Ensembl databases of MANE transcripts (one hg19, the other GRCh38). Without specifying the input, the
database will contain all transcripts in that Ensembl release.

## Option 2b. Make an RefSeq Database

To make an RefSeq database is similar to the Ensembl process.

```
 python3 RefSeqDB.py -h
Usage:

RefSeqDB.py <options>

CAVA refseq v1.3.0

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input filename (list of NM IDs)
  -o OUTPUT, --output=OUTPUT
                        Output filename prefix
  -n, --nm-only         Only collect NM transcripts
  -D OUTPUT_DIR, --outdir=OUTPUT_DIR
                        Output directory
  -r REFSEQ, --release=REFSEQ
                        RefSeq release version
  --no_hg19             Set this to skip hg19 builds

Example usage: CAVA/RefSeq.py -e GCF_000001405.39_GRCh38.p13 -o refseq_db_75
Note: by default, hg19 will be created using crossmap
Version: 1.3.0
```

The big difference here is that you can restrict to only NM (human curated)
transcripts. Note that the RefSeq version system is not as simple as Ensembl. You can find out what version number by
identifying the pattern that fits this expression:
> 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/'
+ version + '/' + version + '_genomic.gtf.gz'

For example, one might want to run the following:

```
python3 RefSeqDB.py -r GCF_000001405.39_GRCh38.p13 -o RefSeq -i data/RefSeq.txt
```
