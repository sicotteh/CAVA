 

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

To install and run CAVA you need:
- Python 3.9 or newer
- GCC and GNU make

Use a virtual environment so CAVA and its Python dependencies stay isolated from your system Python.

4 INSTALLATION ON LINUX OR MAC
------------------------------

### Prerequisites

- Python 3.9+
- `git`
- `git-lfs` (required for downloading published catalog payloads)

Install and initialize Git LFS once per machine:

```bash
git lfs install
```

### Option A: Pull source from GitHub directly (with LFS catalogs)

```bash
git clone https://github.com/sicotteh/CAVA.git
cd CAVA
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Download only the catalog payload you want from Git LFS (example: GENCODE 50 GRCh38):

```bash
git lfs pull --include="catalogs/files/CAVA_GENCODE_50_GRCh38.gz,catalogs/files/CAVA_GENCODE_50_GRCh38.gz.tbi,catalogs/files/CAVA_GENCODE_50_GRCh38.txt,catalogs/files/CAVA_GENCODE_50_GRCh38.cesis"
```

Create a config from the template and set paths:

```bash
cp config_template.txt config.gencode.txt
# Edit config.gencode.txt:
#   @reference = /absolute/path/to/your/reference.fa
#   @ensembl = /absolute/path/to/CAVA/catalogs/files/CAVA_GENCODE_50_GRCh38.gz
```

### Option B: Install from pip/PyPI and fetch catalogs with `cava_data`

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install cava
```

List published catalogs:

```bash
cava_data list
```

Install one catalog payload set into a local data directory:

```bash
cava_data install gencode-50-grch38 ./cava_data --force
```

Create a config file that points `@ensembl` to the installed catalog:

```bash
cava_data config gencode-50-grch38 ./cava_data --output config.gencode.txt --force
# Then edit @reference in config.gencode.txt to your local FASTA path.
```

In a new terminal session, reactivate the same environment with:

```bash
source .venv/bin/activate
```

Windows PowerShell activation:

```powershell
.venv\Scripts\Activate.ps1
```

### Editable installation from source

```bash
git clone https://github.com/sicotteh/CAVA.git
cd CAVA
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Install test dependencies when needed:

```bash
python -m pip install -e ".[tests]"
```

If a source build of `pycurl` is required on your operating system, install the required libcurl development headers first, then rerun `python -m pip install -e .`.


5 RUNNING CAVA
--------------

Before using CAVA, create a config file that provides:

1) A reference FASTA (`@reference`) and matching `.fai`.
2) A transcript catalog (`@ensembl`) from either:
   - GitHub + Git LFS payloads in `catalogs/files`, or
   - `cava_data install ...` output directory.

If you need a quick reference FASTA example:

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
samtools faidx hg38.fa
```

Run CAVA with:

```bash
python -m cava.CAVA -c config.gencode.txt -i input.vcf -o output
```

It requires three command line arguments: 
the name of the configuration file (-c), the name of the input file (-i) 
and the prefix of the output file name (-o). 

The supported Python version for this release is Python 3.9+.

### Phased Haplotype Mode (experimental)

This release adds an experimental mode for phased cis haplotypes encoded in the VCF `ID` field. For example if the ID field is `chr17_7675155_G_A;chr17_7675157_G_C`, then the row-level VCF allele is treated as the contiguous haplotype produced by those atomic edits.

CLI flags:

- `--parseHaplotype` enables parsing/validation of semicolon-separated atomic IDs.
- `--splitBasedOnProtein` emits additional subset records after haplotype annotation when the protein consequence can be partitioned into independent components.
- `--splitAdjacentProtein` emits optional additional records even when the full protein consequence is an adjacent multi-residue replacement, as long as the split records still reproduce the same residue changes after reannotation.

Dependencies:

- `--splitBasedOnProtein` requires `--parseHaplotype`.
- `--splitadjacentprotein` requires `--parseHaplotype`.

Atomic ID encoding:

- Row-level `CHROM`, `POS`, `REF`, `ALT` remain scalar contiguous VCF alleles.
- Row-level `ID` contains atomic components separated by `;`.
- Each atomic token must be `CHROM_POS_REF_ALT` on the genomic plus strand.

Example:

```text
CHROM=chr17
POS=7675155
ID=chr17_7675155_G_A;chr17_7675157_G_C
REF=GCG
ALT=ACC
```

HGVS behavior:

- Adjacent in-phase variants are reported as a merged `delins` at genomic and cDNA HGVS levels.
- Non-adjacent in-phase variants are reported in HGVS allele cis-notation in `HGVSg` and `HGVSc`, with the sequence identifier outside the brackets and the component variants inside the brackets separated by semicolons.
- The `CSN` field remains a single merged cDNA consequence even when `HGVSg` and `HGVSc` use cis-allele notation.
- `HGVSp` continues to reflect the fully reannotated merged protein consequence unless extra split records are emitted.

Split behavior:

- Canonical parsed-haplotype output is emitted when the merged consequence is valid and does not need forced decomposition.
- Split records are identified by the existing provenance tags: `CAVA_ORIGHAPLOTYPE` is the full atomic list and `CAVA_HAPLOTYPE` is the atomic subset used for that emitted line.
- Canonical records have `CAVA_ORIGHAPLOTYPE == CAVA_HAPLOTYPE`.
- Split records have `CAVA_ORIGHAPLOTYPE != CAVA_HAPLOTYPE`.
- Variants may be forcibly split, even without `--splitBasedOnProtein`, when the merged haplotype spans distinct functional regions that should not be represented as one DNA event for downstream interpretation.

Forced-split edge cases:

- If the region between phased variants crosses different functional contexts such as coding sequence, splice region, or non-splice intron, CAVA emits split records instead of a single canonical merged record.
- Variants entirely in UTR sequence that are 1 bp apart are split.
- Adjacent variants entirely in UTR sequence are not forced to split by that rule.
- Existing splice-driven decomposition still applies when the merged splice annotation is not supported by singleton components.

New INFO tags:

- `CAVA_ORIGHAPLOTYPE`: original full atomic token list.
- `CAVA_HAPLOTYPE`: atomic token list used for the emitted record.

Tab-delimited output:

- `HGVSG`, `HGVSC`, and `HGVSP` are followed by `CAVA_ORIGHAPLOTYPE` and `CAVA_HAPLOTYPE` columns.
- The new haplotype columns mirror the INFO-tag content and are present in TSV output for both canonical and split records.

INFO field encoding:

- Semicolons are encoded as `%3B`.
- Percent signs are encoded as `%25`.

Notes:

- Canonical output is usually emitted first, but may be suppressed when forced decomposition is required by splice or mixed-region rules.
- Split outputs are additional records unless canonical output is suppressed.
- If haplotype parsing is requested, malformed IDs or incompatible rows are reported as errors.

6 LICENCE
---------

CAVA is released under MIT licence (see the LICENCE file).

7 CHANGES HISTORY
---------
Version 2.0.15 fixes a few bugs, upgrades packages to newer packages (previous versions were circa 2018), expands unit test coverage, and adds multi-variant haplotype annotation with optional split output when a haplotype can be represented as independent consequences.
This version of CAVA includes the following changes (aside from bug fixes, especially for edge cases where multiple interpretations could apply)
- support refseq transcripts in addition to ensembl
- include new tags: CAVA_HGVSG, CAVA_HGVSC, CAVA_HGVSP to represent the current full HGVS nomenclature (G=Genomics, C=CDNA,P=Protein) for the HGVSC and HGVSP. We do not support the genomic tandem repeats for HGVSG nor imperfect repeats. These fields must be URL-decode (uudecode) because they include ';' encoded as %3B (';' is not a legal character in the VCF INFO field).
   -- This requires the addition of a files to support the protein information. Catalogs prior to version 2.03 will not be compatible because of that additional file.
- parsed haplotypes now follow HGVS more closely:
   -- adjacent phased variants use merged `delins` at genomic and cDNA HGVS levels
   -- non-adjacent phased variants use allele cis-notation in `HGVSg` and `HGVSc`
   -- `CSN` remains a merged cDNA consequence for parsed haplotypes
- parsed haplotypes now decompose automatically when singleton components fall into different functional regions, including coding versus splice/intron and UTR variants separated by 1 bp
- Include support for selenocysteine genes and alternate stop codon. Usually this is a conservative interpretation, often resulting in p.? when a novel stop codon might be discovered.
- Support for MANE transcripts (including transcripts with UTR's of length 0  which are part of MANE 1.0)
- Includes expanded unit test coverage for haplotype HGVS formatting, split provenance, splice/coding mixed decomposition, and UTR spacing edge cases.
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



