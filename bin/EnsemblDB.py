#!env/bin/python3

from optparse import OptionParser
from ensembldb import main
import os
with open(os.path.join('bin', 'VERSION')) as version_file:
    version = version_file.read().strip()


# Command line argument parsing
descr = 'bin/EnsemblDB.py' + version
epilog = '\nExample usage: bin/EnsemblDB.py -e 75 -o ensembl_db_75\n' \
        'Note: by default, hg19 will be created using crossmap\n' \
        'Version: {}\n' \
         '\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nbin/EnsemblDB.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--input", default=None, dest='input', action='store', help="Input filename (list of ENST IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
parser.add_option('-e', "--ensembl", default=None, dest='ensembl', action='store', help="Ensembl release version.py")
parser.add_option("--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")

(options, args) = parser.parse_args()

if int(options.ensembl) <= 75:
    options.no_hg19 = False
    print('Since this is <= version 75, no hg19 conversion will be done.')

options.select = False
if options.input is not None:
    options.select = True

options.version = version

main.run(options)