# Preparing a transcript database

These instructions have two parts. Steps 1, 2's, and 3 are the basic instructions. 
Following that, there are two special sections. One for hg19/GRCh37 catalogs and one for 
non-human catalogs. Note that by default, the instructions in section 2 can automatically create
an hg19 catalog. These 'automatic' hg19 catalogs are based on using cross-map for mapping transcripts. 
Catalogs created using cross-map may have issues if the reference sequence of the exons is not sufficiently
identical between the two builds to avoid issue (no gaps and presence of start and stop codon, without additional
in-frame stop codons). The instructions in section 4 avoids those issues.

## Step 1. Identify a transcript set (optional for Refseq, required for Ensembl)

The best place to get a list of transcript to use is from the [MANE](ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/)
project. These are curated & versioned lists of what should be considered the "default" transcript for most genes.

In this case, we will use release 1.3, but other options are listed further down.

```
# Download GTF files for MANE 
export MANEVERSION="1.4"
export ENSEMBLVERSION="114"

wget -O data/MANE${MANEVERSION}_ENST.gtf.gz https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANEVERSION}/MANE.GRCh38.v${MANEVERSION}.ensembl_genomic.gtf.gz 
wget -O data/MANE${MANEVERSION}_RefSeq.gtf.gz https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANEVERSION}/MANE.GRCh38.v${MANEVERSION}.refseq_genomic.gtf.gz
# Download GTF files for current Refseq.
wget -O data/RefseqLatest_GRCh38_Refseq.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
# Download GTF files for latest GRCh37 annotation release
wget -O data/Refseq_GRCh37_Refseq.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz
# Download GTF file for build 38 ensembl release
wget -O data/Homo_sapiens.GRCh38.114.gtf.gz https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
# Download GTF file for build 37 ensembl release
wget -O data/Homo_sapiens.GRCh37.87.gtf.gz https://ftp.ensembl.org/pub/grch37/release-${ENSEMBLVERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz
wget -O data/Homo_sapiens.GRCh37.75.gtf.gz https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip -c data/Homo_sapiens.GRCh37.75.gtf.gz | wc -l
#2828317
gunzip -c data/Homo_sapiens.GRCh37.87.gtf.gz | wc -l
#2612766

gunzip -c data/Homo_sapiens.GRCh37.75.gtf.gz | wc -l
#2828317
gunzip -c data/Homo_sapiens.GRCh37.87.gtf.gz | wc -l
#2612766

```


### Create lists of wanted transcripts
```
# MANE Ensembl transcripts.
zcat data/MANE${MANEVERSION}_ENST.gtf.gz |gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}' | perl -ane 'if ($_ =~ /transcript_id "([^"]+)"/){print $1 . "\n"}' |grep ENST|sed 's/;//;s/\"//g'|sort -u > data/MANE${MANEVERSION}_ENST.txt
#All Recent build 38 transcripts.
zcat data/Homo_sapiens.GRCh38.${ENSEMBLVERSION}.gtf.gz  |gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}' | perl -ane 'if ($_ =~ /transcript_id "([^"]+)"/){print $1 . "\n"}' |grep ENST|sed 's/;//;s/\"//g'|sort -u > data/Ensembl${ENSEMBLVERSION}.GRCh38.ENST.txt
zcat data/Homo_sapiens.GRCh37.87.gtf.gz  |gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}' | perl -ane 'if ($_ =~ /transcript_id "([^"]+)"/){print $1 . "\n"}' |grep ENST|sed 's/;//;s/\"//g'|sort -u > data/Ensembl87.GRCh37.ENST.txt
zcat data/Homo_sapiens.GRCh37.75.gtf.gz  |gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}' | perl -ane 'if ($_ =~ /transcript_id "([^"]+)"/){print $1 . "\n"}' |grep ENST|sed 's/;//;s/\"//g'|sort -u > data/Ensembl75.GRCh37.ENST.txt
#Refseq transcripts, limited to NM_ (e.g. ) with a protein..

zcat data/MANE${MANEVERSION}_RefSeq.gtf.gz | gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}'| perl -ane 'if ($_ =~ /transcript_id "([^"]+)"/){$NM=$1;print $NM . "\n"}' | uniq | sort | uniq > data/MANE_${MANEVERSION}.RefSeq.txt
```
Some transcripts are mapped on alternate contigs instead of the reference genome. For example, 
Refseq transcript NM_014219.3 is NOT derived from the reference genome .. and the genomic sequence does
not have a proper translated protein where this transcripts maps to. 





The MANE GTF also defines transcripts on alternate 
sequences .. and they give them alternate names.  e.g. NM_014219.3_3 is annotated and produces a protein. 
Since CAVA generate transcripts from the underlying sequence, these alt-mapped transcripts may not produce the same mRNA as the base reference transcript.
!!! Including those may produce errors if a sequence-based stop codon is encountered before the expected end of the annotated transcript!!!
Ideally, exclude from the catalog by limiting to GRCh38 1-22,X,Y,M

Reviewing transcripts with an exception tag, we found that transcript NM_001424184.1 is missing one base to actually get a stop codon, and the next base in the genome is not an "A". This A is added by polyadenylation and creates a stop codon. It has to be excluded as it is partial.

Annotation for 'alternative start codon' do not matter as long as the GTF properly points to it. This version 2.0.14 of CAVA or later only.
AUU and AUA are accepted alternative start codons (Normal Methionine is AUG .. so only the wobble position changes).


For a few genes, the reference sequence sequence does not include the gene, they are instead located on alt or fixed contigs.
This command will show which genes (those without an '_' in their name) are primarily mapped on an alt contig.

gunzip -c data/MANE${MANEVERSION}_RefSeq.gtf.gz | gawk -F "\t" '$3=="transcript"{print $0}' | perl -aE '@a=split(/\t/,$_);if ($a[0] =~ /_/){if($_ =~ /gene_id "([^"]+)"/) {if(! ($1 =~ /_/)) {$gene=$1;$_ = ~ /transcript_id "([^"]+)"/; $transcript = $1;print $a[0] . "\t" . $gene . "\t" . $transcript . "\n"}}}' |  uniq

chr1_KQ031383v1_fix	PRAMEF22	NM_001100631.2
chr1_MU273333v1_fix	LOC102724250	NM_001405530.2
chr11_JH159136v1_alt	OR8U8	NM_001013356.2
chr11_JH159136v1_alt	OR8U9	NM_001013357.1
chr11_JH159137v1_alt	OR9G9	NM_001013358.2
chr12_GL877876v1_alt	TAS2R45	NM_176886.2
chr19_GL949746v1_alt	LILRA3	NM_006865.5
chr22_KI270879v1_alt	GSTT1	NM_000853.4
chr17_KI270909v1_alt	CCL3L1	NM_021006.6
chr17_KI270909v1_alt	CCL4L1	NM_207007.4
chr6_GL000254v2_alt	LOC110384692	NM_001352000.1


The GSTT1 can be added to the catalog. LILRA3,CCL3L1,CCL4L4 also have variants in Clinvar. The other genes are not as important.



The Ensembl Catalog has a lot more transcripts that have multiple location (both ref and alt contigs)
The only transcript we have to exclude is : ENST00000434431.2 because it is partial (this is NM_001424184.1 , described above)

gunzip -c data/MANE1.4_ENST.gtf.gz | gawk -F "\t" '$3=="transcript"{print $0}' | perl -aE '@a=split(/\t/,$_);if ($a[0] =~ /_/){if($_ =~ /gene_name "([^"]+)"/) {if(! ($1 =~ /_/)) {$gene=$1;$_ = ~ /transcript_id "([^"]+)"/; $transcript = $1;print $a[0] . "\t" . $gene . "\t" . $transcript . "\n"}}}' |  uniq | cut -f2 > ensembls.alt

for i in `cat ensembls.alt `; do  n=`gunzip -c data/MANE1.4_ENST.gtf.gz | grep $i | gawk -F "\t" '$3=="transcript"' | grep -v '_alt' | grep -v '_fix' | wc -l`; if [[ $n -eq 0 ]] ; then echo "$i";fi;done


# Skip all transcripts with pseudo tag
gunzip -c data/RefseqLatest_GRCh38_Refseq.gtf.gz | gawk -F "\t" '$3=="start_codon"{print $0}' | perl -aE '@a=split(/\t/);if($a[8] =~ /gene_id "([^"]+)"/) {$gene=$1;$a[8] =~ /transcript_id "([^"]+)"/; $transcript = $1;if ($transcript =~ /_.*_/){print $_}}'  | grep -v 'unassigned_transcript' | grep 'pseudo "true' | more



The ABO gene is not part of MANE, but, given it's importance, you may want to add it. It is not known to produce a protein, but it is transcribed. THat is the only
requirement for safely adding a gene.
Only 1 refseq transcript is defined for ABO
NM_020469.3  which matches ENST00000644422.3
```
echo "NM_020469.3" >> data/MANE${MANEVERSION}_RefSeq.txt
echo "ENST00000644422.3" >> data/MANE_${MANEVERSION}_ENST.txt

```


## Step 2
You have to create a transcript database. You have 4 choices.
- Ensembl
- Refseq,
- ensembl, or
- MANE (which limits to clinically relavant transcripts .. mostly 1 transcript/gene, but ~ 100 genes have multiple transcripts
    MANE collection can be either RefSeq transcript ids or Ensembl transceript ids   
    If you choose the MANE subset .. you might also want to create either the "All ensembl" or the "All Refseq" in order to add some genes missing from MANE (in particular ABO)

### Step 2a. Make an Ensembl Database

To make an Ensembl database, use the following program:

```
python3 $CAVAGITREPO/CAVA/cava/EnsemblDB.py -h
```

``` Usage: 

EnsemblDB.py <options>

EnsemblDB.py2.0.13

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
                        Ensembl release version.py (build GRCh37 for
                        build<=75)/
  -x, --no_hg19         Set this to skip hg19 builds
  -b BUILD, --build=BUILD
                        GRCh38 or something else. If GRCh38, then will get
                        remapped to hg19 unlesss --no_hg19 isset
  -u URL_GTF, --url_gtf=URL_GTF
                        Download an arbitrary GTF. Overrides ensembl option.
                        e.g. use to get build 111 GRCh37: https://ftp.ensembl.
                        org/pub/grch37/release-
                        111/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz

Example usage: EnsemblDB.py -e 75 -o ensembl_db_75
Note: by default, hg19 will be created using crossmap
Version: 2.0.13
```


By default, if you specify a version > 75, the coordinates will be in GRCh38. Additionally, this script will
automatically create an hg19 coordinate version unless you tell it not to (--no_hg19 option).
Using our example from above, the command

```
python3 $CAVAGITREPO/CAVA/cava/EnsemblDB.py -e ${ENSEMBLVERSION} -o ENST101 -D data -i data/ENST.txt
``` 


will create two Ensembl databases of all the transcripts in data/ENST.txt (which were the MANE transcripts).
One GRCh38 and one hg19 transcripts databases.
Without specifying the input (-i) , the database will contain all transcripts in that Ensembl release, including non-human ones.

### Step 2b. Make an RefSeq Database

To make an RefSeq database is similar to the Ensembl process.

```
python3 RefSeqDB.py -h
Usage: 

RefSeqDB.py <options>

CAVA refseq v2.0.13

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
                        This is no longer used to create url for download.
                        Browse https://ftp.ncbi.nih.gov/genomes/refseq/vertebr
                        ate_mammalian/Homo_sapiens/annotation_releases/  or ht
                        tps://ftp.ncbi.nlm.nih.gov//genomes/refseq/vertebrate_
                        mammalian/Homo_sapiens/all_assembly_versions
  -x, --no_hg19         By default remap to hg19. Set this option if
                        downloadin GRCh37 directly or for non-human organisms
  -b BUILD, --build=BUILD
                        GRCh38 or something else. If GRCh38, then will get
                        remapped to hg19 unlesss --no_hg19 isset
  -u URL_GTF, --url_gtf=URL_GTF
                        Download an arbitrary GTF.

Example usage: RefSeqDB.py  -u 'https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz' -o refseq_db_110_GRCh38 -D refseq_110
Note: by default, hg19 will be created using crossmap
Version: 2.0.13

```

For CAVA, you are encouraged to restrict to only NM transcripts. 
NM transcripts are always accurately transcribed from the reference genome and have a complete CDS. 
Predited transcripts can cause issues with CAVA.



### Step 2c - Make a MANE database

```
python3 MANE.py -h
```
```
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -D OUTPUT_DIR, --outdir=OUTPUT_DIR
                        Output directory
  -e version, --mane_version=version
  --no_hg19             Set this to skip hg19 builds

Example usage: MANE.py -e ${MANEVERSION} -o mane_${MANEVERSION}
Note: by default, hg19 will be created using crossmap.
Version: 1.3.3
```


You should add the ABO transcript if you want it (extracting it from the refseq catalog)



## Step 3
Finally, we have to create annotation for selenocysteine genes
```
cd $CAVAGITREPO/CAVA/cava/ensembldb
```

If you are only interested in the latest refseq/MANE data on build GRCh38, you can skip the Step 3a
If you are interested in either the ensembl data or the GRCh37 data, the following may be needed.
### Step 3a
In Step 3b, we will search for transcripts with CECIS annotation in the NCBI Entrez database.
The NCBI entrez database only allows to search for annotation on the latest refseq version.
However once we know that a refseq has CESIS annotation, we can fetch earlier refseq versions
corresponding to older builds or alternate matches.

To get a mapping file between ensembl and refseq in GRCh37, use the last build of GRCh37 in ensembl.
In fact, we already created this file in the distribution. Two columns,first one is an ensembl transcript ID, the second one is the full refseq.version

```
mysql -u anonymous  -h ensembldb.ensembl.org  -D homo_sapiens_core_75_37 -e 'SELECT transcript.stable_id, xref.display_label FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = "Transcript"  AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = "RefSeq_mRNA"' > ensembl_refseq_GRCh37_75.txt
```

If you want the latest GRCh38 ensembl, you can do the following. At this time, releast ${ENSEMBLVERSION} is the latest. Update the '${ENSEMBLVERSION}' as needed.

```
mysql -u anonymous  -h ensembldb.ensembl.org  -D homo_sapiens_core_${ENSEMBLVERSION}_38 -e 'SELECT transcript.stable_id, xref.display_label FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = "Transcript"  AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = "RefSeq_mRNA"' > ensembl_refseq_GRCh38_${ENSEMBLVERSION}.txt
```

the first few lines are

stable_id	display_label
ENST00000612465	NM_001347681.2




Note that different ensembl can map to the same refseq in different build. Multiple ensembl transcripts can map to the same refseq in the same build. 

### Step 3b
The following command queries NCBI for selenocysteine proteins and their annotation.
By specifying the --ensembl_refseq option, it allows the annotation of ensembl 
transcripts by getting their data from the equivalent refseq transcript.

```
python3 create_selenodata.py  -h
```

```
Usage: 
Example usage: Run after catalog file were created: Without any option, will build for latest Refseqs for Human.
 python3 create_selenodata.py 


Options:
  -h, --help            show this help message and exit
  -e ENSEMBL_REFSEQ, --ensembl_refseq=ENSEMBL_REFSEQ
                        tab-delimited file of ensembl and matching
                        refseqs.versiom. If specified, overrides the --refseqs
                        option
  -r REFSEQS, --refseqs=REFSEQS
                        file with refseq.version, the version number can make
                        the refseqs build-=specific.
  -D OUTPUT_DIR, --outdir=OUTPUT_DIR
                        Output directory
  -g GENES, --genes=GENES
                        file containes list of known selenocysteine genes
                        (organism-specific)
  -o ORGANISM, --organism=ORGANISM
                        Organism to use for the entrez query (defaults to
                        homo_sapiens).  No spaces in name (e.g.
                        mus_musculus,sus_scrofa)
  -t BUILD_TAG, --tag=BUILD_TAG
                        arbitrary tag to put in output filenames. Defaults  to
                        refseq or ensembl.

```

Build selenocysteine database. The list of NM_ numbers (with the latest version number) with SECIS annotation is discovered
by querying entrez. By providing a --refseqs option, it allows earlier version numbers to
be annotated. The following will create a file data/SECIS_in_refseq_pos.homo_sapiens.refseq_GRCh38_latest.txt
```
python3 create_selenodata.py  --refseqs data/RefSeq.txt -O data -t refseq_latest
```

Create Selenoysterine data, but replace refseq ids by Entrez Stable ids (without a .version). Entrez will be queried for the
refseq version specified in the 2nd column of the file (if a refseq with different/same version is known to have SECIS annotation). <br>
The following will create a file data/SECIS_in_refseq_pos.homo_sapiens.ensembl_GRCh38_${ENSEMBLVERSION}.txt

```
python3 create_selenodata.py  --ensembl_refseq ensembl_refseq_GRCh38_${ENSEMBLVERSION}.txt -O data -t ensembl_GRCh38_${ENSEMBLVERSION}
```
The following will create a file data/SECIS_in_refseq_pos.homo_sapiens.ensembl_GRCh37_75.txt

```
python3 create_selenodata.py  --ensembl_refseq ensembl_refseq_GRCh37_75.txt -O data -t ensembl_GRCh37_75
```

##### This code will generate a config file to properly annotate selenocysteine genes:
<B>Biological mechanism:</B> UGA codons (normally a  codon) ealier than 51 bp 5' from the SECIS element (usually located in the 3' UTR, but possibly in the coding region)
can be recoded as Sel. 
THere is some uncertainty as to the range of action of the SECIS element. 
The range of activity could be as far as 111 bases.

The config file will contain the location of the 3'most SECIS element 
(some genes have multiple) as well as the location of the latest (3'-most) Selenocysteine 
(UGA recoded as Sel). The position of the known Sel is a position that guarantees that 
any UGA before that position will be recoded as a Sel.


This code queries NCBI for selenocysteine transcripts by using known genes as well as genes with SECIS elements.
It creates the following file.
SECIS_in_refseq_pos.homo_sapiens.txt  (the words homo_sapiens will be whatever your organism is)
If you specify the --ensembl_refseq option, it will also create an ensembl file.
SECIS_in_ensembl_pos.homo_sapiens.txt  

You can get appropriate files for different build by specifying --ensembl_refseq or -refseqs files. 
The refseq version number (number after the period in a refseq.. e.g. the '4' NM_001127670.4 
is the latest GRCh38 for KCNE1 gene and NM_001127670.1 (version '1') is the most recent 
version on GRCh37. This paper(PMID=34129815), says that 790 transcripts with nucleotide differences between
those two builds (transcripts not all reported), so most transcripts should have identical coding region
between builds. The transcript versions differences tend to affect the UTR length, so our code
doesn't use absolute transcript coordinates to deal with SECIS elements, but keeps track of them 
relative to the start codon.

By default, the program will scan for any and all genes that appear to have SECIS 
annotation in addition to using a hardcoded gene list for humans.
If you specify a different organism, it is strongly suggested that you do a little research and find out the gene symbols for selenocysteine genes in your organism
and add the -g gene_list option.

Note that can only specify one of the --ensembl or the --refseqs option.



Copy this file along with your catalog to a file named (whatever your catalog name)

```
cp CECIS_in_refseq_pos.txt CATALOGNAME.cesis
```

The file will contain the following fields: (shown here with refseq, but would be ensembl with  --ensembl_refseq option)
gene	SECIS_maxstart	SECIS_maxend	accn	last_sec_pos	cds_start	cds_end	mRNA_len
SELENOP	1497	1579	NM_001085486.3	1231	100	1245	2085
SELENOP	1556	1638	NM_001093726.3	1290	69	1304	2144

### Build 37 Config Files

#### Refseq Build 37
Browse this <br \>
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions

The latest annotation for build 37 would be gotten with a command like<br />
```
python3 RefSeqDB.py --no_hg19 -b GRCh37 --nm-only --url_gtf "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
```

#### MANE transcripts for Build 37
The MANE transcripts are defined based on GRCh38. We can still get the GRCh37 transcripts, but
they will be different version numbers and so might have diffferent UTR's. Refseqs sequence are defined to match
the nucleotide sequence of the build they are on. So there may be nucleotide difference of a build-38 refseq version
compared to a build 37 version. In the few cases where these differences prevent the protein from being assembled, the refseq protein
ids are set to the empty string in the catalog annotation.(1 known case)
We cannot guarantee that all MANE transcripts will be available, but we are using MANE to pick the subset of clinically relevant isoforms.

There are 4 steps. 
1. Get the build 38 refsseq Transcripts ids and specific build 38 version.
2. Get the available build 37 transcripts and version from the most recent GTF
3. Create a file of the build 37 transcripts version corresponding to a MANE transcript.
4. Pick out the subset of the step 2 GTF matching the MANE transcripts.

The only complicated this is that some Refseq transcripts have a 'madeup' id with an extra '_',
such as NM_000595.4_6. This indicates that that the transcript sequence is from build 38 and
that the sequence is different than the build 37 sequence. When the sequence does not match, 
there may be an entry with the unmodified refseq id , but the protein id will be empty.

#### Steps 1-3
```
wget -O data/MANE.GRCh38.v${MANEVERSION}.refseq_genomic.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANEVERSION}/MANE.GRCh38.v${MANEVERSION}.refseq_genomic.gtf.gz
wget -O data/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz

# All GRCH37 transcipts for which we have proper annotation on chrs1-22,X,Y,M
gunzip -c data/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz | gawk -F "\t" '$1 ~ /^NC_[^_]+$/{print $9}' | grep transcript | perl -ane 'if($_ =~ /transcript_id "([^"]+)"/) {print $1 . "\n"}' | uniq | sort | uniq > Refseq_GRCh37_transcripts.txt
# All MANE GRCh38 transcripts, irrespective of wether they map properly (because transcripts with sequence differences on GRCh38 may match the GRCh37 sequence)
gunzip -c data/MANE.GRCh38.v${MANEVERSION}.refseq_genomic.gtf.gz | gawk -F "\t" '$1 ~ /^[chr]*[^_]+$/{print $9}' | grep transcript | perl -ane 'if($_ =~ /transcript_id "([^"]+)"/) {print $1 . "\n"}' | uniq | sort | uniq > MANE${MANEVERSION}.GRCh38.transcripts.txt
gawk -F "\t" 'BEGIN{while(getline<"MANE${MANEVERSION}.GRCh38.transcripts.txt"){MANEV[$1]=1;split($1,a,".");MANE[a[1]]=1}}MANEV[$1]==1{print $1;next}{split($1,a,".");if(MANE[a[1]]==1){print $1}}' Refseq_GRCh37_transcripts.txt | sort | uniq | perl -ane 'if(! ( $_ =~ /\.[0-9]_/)){print $_}' > tokeep.GRCh37
```

Watch the error messages. We can see that one of the transcript of MANE (NM_004892.6) 
does not end up in the catalog since there is no protein_id in the gtf because of significant frameshifts between 
the transcript sequence (from GRh38) and GRCh37.


In MANE ${MANEVERSION}, some of the refseq (also) have the ids that are not purely refseq, they are derived from refseq 
ids with an extra  '_[0-9]' suffix. For MANE 1.3 this is what we got.

- 19290 transcripts in GRCh38 MANE 1.3 on chroms 1-22,X,Y,M
- 19171 matching refseq ids with build 37 annotation (file tokeep.GRCh37), including some GRCh38 refseqs with a suffix if none was available without the suffix


#### Step 4, final catalog creation.
```
cd ../..
python3 cava/RefSeqDB.py --no_hg19 -b GRCh37 --nm-only -i cava/ensembldb/tokeep.GRCh37  -D cava/ensembldb/data -o GRCh37.MANE${MANEVERSION}.refseq  --url_gtf "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
```

### GRCh37 MANE with ensembl transcripts with newest annotation from refseq
To get the build GRCh37 MANE ensembl transcripts, create the refseq file (previous section)
and just use the ensembl_refseq_GRCh37.txt to substitute the ids

We have to allow unexact match (e.g. different version), otherwise we only get 491 transcripts
```
# 
gunzip -c cava/ensembldb/data/GRCh37.MANE${MANEVERSION}.refseq.gz | gawk -F "\t" 'BEGIN{while(getline<"cava/ensembldb/ensembl_refseq_GRCh37_75.txt"){if ( $1 ~ /ENST.*/) {enst=$1;refseq=$2;split(refseq,arr,".");refseq0=arr[1];if(refseq in R2E){split(R2E[refseq],arr,",");arr[length(arr)+1]=enst;elem=arr[1];for(i=2;i<=length(arr);i+=1){elem = elem "," arr[i]};R2E[refseq]=elem} else{R2E[refseq]=enst}; if(refseq0 in R02E){split(R02E[refseq0],arr,",");arr[length(arr)+1]=enst;elem=arr[1];for(i=2;i<=length(arr);i+=1){elem = elem "," arr[i]};R02E[refseq0]=elem} else{R02E[refseq0]=enst};  }}} $0 ~ /#.*/{print $0;next}{transcript = $1;split($1,arr,".");nm=arr[1];if(transcript in R2E){n=split(R2E[transcript],arr,",");for(i=1;i<=n;i+=1) {OFS="\t";$1=arr[i];print $0}} else {if(nm in R02E){n=split(R02E[nm],arr,",");for(i=1;i<=n;i+=1) {OFS="\t";$1=arr[i];print $0}} }}' > cava/ensembldb/data/MANE_ensembl_GRCh37.gtf
bgzip cava/ensembldb/data/MANE_ensembl_GRCh37.gtf
```


### GRCh37 ***MANE*** with old annotation. (can use ensembl 75 .. or provide the url and download 87/111 annotaion)
```
cd $CAVAGITREPO/CAVA
python3 cava/EnsemblDB.py -e 75 -o ENST75_GRCh37 -D ensembldb/cava/data --no_hg19 -b GRCh37 -i cava/ensembldb/data/MANE${MANEVERSION}_ENST.txt
``` 
or
```
cd $CAVAGITREPO/CAVA/
python3 cava/EnsemblDB.py -o ENST87_GRCh37 -D data  -e 87 -b GRCh37 --no_hg19 -i cava/ensembldb/data/MANE${MANEVERSION}_ENST.txt --url_grf "https://ftp.ensembl.org/pub/grch37/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
``` 
### GRCh37 ***ensembl*** with old annotation. (can use ensembl 75 .. or provide the url and download 87/111 annotation)
```
cd $CAVAGITREPO/CAVA
python3 cava/EnsemblDB.py -e 75 -o ENST75_GRCh37 -D ensembldb/cava/data --no_hg19 -b GRCh37 
``` 
or
```
cd $CAVAGITREPO/CAVA/
python3 cava/EnsemblDB.py -o ENST87_GRCh37 -D data  -e 87 -b GRCh37 --no_hg19 --url_grf "https://ftp.ensembl.org/pub/grch37/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
``` 



### GRCh37 ***ensembl*** annotation for all the transcripts included in the most recent Ensembl Build 38 release.
An ensembl GRCh37 catalog will be the least up-to-date because the the ensembl transcript annotation for GRCh37 is frozen. 
Even if you can find a more recent annotation file (upgraded reference), the gene/transcript annotation is old.
However, you can download the refseq annotation and swap the refseq transcript id for it's matching ensembl ids to get updated annotation.

Even though some ensembl annotation is available listed as build 111 or build 87, the transcript annotation/mapping was really made on built 75.
The following approach will get the transcripts annotation for a current ensembl release.
- Identify wanted ensembl release and define R38 variable 
    - export R38=111
- Create a build 38 CAVA catalog (with mapping) from a specific catalog, and allow build 37 mapping.
    - Create a list of the recent build 38 ensembl transcripts ids --> Ensembl_GRCh38.${RELEASE}.ids (LIST_E38_${R38})
- Find build-37-native transcript annotation from refseq
    - Create a Refseq catalog for build 37 using the latest annotation release. 
    - create a list of transcripts: LIST_R37_LATEST
    - Download the refseq<->ensembl mappings for build 38 (LIST_M38_${RELEASE})
    - Swap the ensembl<->refseq ids in the catalog files 
    - subset catalog files down to overlap between mapping files LIST_R37 and LIST_X38
- Create a build-37-native transcript CAVA catalog from ensembl (release 75/86/111 --> RELEASE37)
    - LIST_E37_${RELEASE37}
    - Extract entries (if present) from LIST_E38_R${REL38} - (LIST_R37^LIST_X38)
- Lastly, for any remaining ids, use the mapping option with the ensembl code.

gunzip -c data/Homo_sapiens.GRCh37.87.gtf.gz | grep  protein_id | grep -v  cds_start_NF | grep -v cds_end_NF | perl -ane 'if($_ =~ /gene_biotype \"protein_coding/){print $_}' |  grep -v polymorphic_pseudogene | grep -v non_stop_decay | grep -v nonsense_mediated_decay | grep -v '#' | gawk -F "\t" '($1=="X"  ||  $1 == "Y" ||  $1 == "MT" || $1 ~ /^[0-9]+$/){print $0}' | perl -ane 'if ($_ =~ /transcript_id \"([^\"]+)\"/){print $1 . "\n"}' |uniq | sort | uniq > 87.core.list

python3 cava/EnsemblDB.py -e 111 -x -o ensembl_37 -b GRCh37 -u 'https://ftp.ensembl.org/pub/grch37/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz'

