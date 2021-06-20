#!env/bin/python3
import os
from optparse import OptionParser
from cava.ensembldb import main_refseq as main

with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

# Command line argument parsing
descr = 'CAVA refseq v' + version
epilog = '\nExample usage: RefSeqDB.py  -e GCF_000001405.39_GRCh38.p13 -o refseq_db_75 -D refseq_75\n' \
        'Note: by default, hg19 will be created using crossmap\n' \
        'Version: {}\n' \
         '\n'.format(version)

OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nRefSeqDB.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--input", default=None, dest='input', action='store', help="Input filename (list of NM IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-n', "--nm-only", default=False, dest='nm_only', action='store_false', help="Only collect NM transcripts")
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
parser.add_option('-r', "--release", default=None, dest='refseq', action='store', help="RefSeq release version")
parser.add_option("--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")

(options, args) = parser.parse_args()

if 'GRCh38' in options.refseq and options.no_hg19 is False:
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
