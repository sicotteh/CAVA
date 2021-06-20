from optparse import OptionParser
from cava.ensembldb import mane_db_prep as main
import os
with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

# Command line argument parsing
descr = 'MANE.py' + version
epilog = '\nExample usage: MANE.py -e 0.91 -D mane_0.91\n' \
        'Note: by default, hg19 will be created using crossmap\n' \
        'Version: {}\n' \
         '\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nMANE.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
parser.add_option('-e', "--mane_version", default=None, dest='ensembl', action='store', help="release version")
parser.add_option("--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")

(options, args) = parser.parse_args()

options.select = False

options.version = version
if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

main.run(options)
