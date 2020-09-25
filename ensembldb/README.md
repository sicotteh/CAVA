# Preparing a transcript database

```
python3 CAVA/bin/EnsemblDB.py -h

Usage: CAVA-1.2.5/ensembl_db <options>

CAVA ensembl_db v1.2.5

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input filename (list of ENST IDs)
  -o OUTPUT, --output=OUTPUT
                        Output filename prefix
  -e ENSEMBL, --ensembl=ENSEMBL
                        Ensembl release version.py
  --no_hg19             Set this to skip hg19 builds (will build by default)
  -s, --select          Select transcript for each gene [default: False]

Example usage: CAVA/ensembl_db -e 75 -o ensembl_db_75
Note: by default, hg19 will be created using crossmap
Version: 1.2.5

```
