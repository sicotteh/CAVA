#!env/bin/python3
import os
import sys
from optparse import OptionParser
from cava.ensembldb import main_refseq as main


with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

# Command line argument parsing
descr = 'CAVA refseq v' + version
epilog = '\nExample usage: RefSeqDB.py  -u ' + \
         '\'https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz\' ' + \
        '-o refseq_db_110_GRCh38 -D refseq_110\n' + \
        'Note: by default, hg19 will be created using crossmap\n' + \
        'Version: {}\n' \
         '\n'.format(version)

OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nRefSeqDB.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--input", default=None, dest='input', action='store', help="Input filename (list of NM IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-n', "--nm-only", default=False, dest='nm_only', action='store_false', help="Only collect NM transcripts")
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")

parser.add_option('-r', "--release", default=None, dest='refseq', action='store',
                  help="This is no longer used to create url for download. Use --url_gtf option\n" +
                  "Browse https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/"
                  +"\n or https://ftp.ncbi.nlm.nih.gov//genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions")
parser.add_option('-x',"--no_hg19",  action='store_false', default=True, dest='no_hg19', help="By default remap to hg19. Set this option if downloadin GRCh37 directly or for non-human organisms")
# latest GRCh37 is here: https://ftp.ncbi.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz

parser.add_option('-b', "--build", default='GRCh38', dest='build', action='store', help="GRCh38 or something else. If GRCh38, then will get remapped to hg19 unlesss --no_hg19 isset")
parser.add_option('-u',"--url_gtf",  action='store', default=None, dest='url_gtf', help="Download an arbitrary GTF.")

(options, args) = parser.parse_args()

if options.refseq is not None:
    sys.stderr.write("--release/-r option is no longer supported as refseq urls are changing. Please use --url_gtf option . " +
                     "\nhttps://ftp.ncbi.nlm.nih.gov//genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions" + " or \n" +
                     "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/\n" )
    exit(1)

if 'GRCh38' in options.build and options.no_hg19 is False:
    print('\nNo hg19 conversion will be done.')

if not options.nm_only:
    print('\nOnly printing out NM transcripts.')

options.select = False

if options.input is not None:
    options.select = True

if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

options.version = version

main.run(options)
