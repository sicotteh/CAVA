#!env/bin/python3

from optparse import OptionParser
from cava.ensembldb import main
import os
import sys
with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

# Command line argument parsing
descr = 'EnsemblDB.py' + version
epilog = '\nExample usage: EnsemblDB.py -e 75 -o ensembl_db_75\n' \
        'Note: by default, hg19 will be created using crossmap\n' \
        'Version: {}\n' \
         '\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nEnsemblDB.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--input", default=None, dest='input', action='store', help="Input filename (list of ENST IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
parser.add_option('-e', "--ensembl", default=None, dest='ensembl', action='store', help="Ensembl release version.py (build GRCh37 for build<=75)/ ")
parser.add_option('-x',"--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")
parser.add_option('-b', "--build", default='GRCh38', dest='build', action='store', help="GRCh38 or something else. If GRCh38, then will get remapped to hg19 unlesss --no_hg19 isset")
parser.add_option('-u','--url_gtf',  action='store', default=None, dest='url_gtf', help="Download an arbitrary GTF. Overrides ensembl option. e.g. use to get build 111 GRCh37: https://ftp.ensembl.org/pub/grch37/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz")


(options, args) = parser.parse_args()
if options.url_gtf is not None:
    if options.ensembl is None:
        sys.err.write("You also need to specify options.ensembl even if using --url_gtf option\n")
        exit(1)
    if options.build is None:
        sys.err.write("You also need to specify options.build  (GRCh38 or GRCh37) even if using --url_gtf option\n")
        exit(1)

if int(options.ensembl) <= 75:
    options.no_hg19 = False
    print('Since this is <= version 75, no hg19 conversion will be done.')

if options.build == "GRCh37":
    options.no_hg19 = True

options.select = False
if options.input is not None:
    options.select = True

options.version = version
if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

main.run(options)