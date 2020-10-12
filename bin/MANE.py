from optparse import OptionParser
from ensembldb import mane_db_prep as main

with open('/CAVA/VERSION') as version_file:
    version = version_file.read().strip()


# Command line argument parsing
descr = 'bin/MANE.py' + version
epilog = '\nExample usage: bin/MANE.py -e 0.91 -o mane_0.91\n' \
        'Note: by default, hg19 will be created using crossmap\n' \
        'Version: {}\n' \
         '\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='\n\nbin/MANE.py <options>', version=version, description=descr,
                      epilog=epilog)
parser.add_option('-D', "--outdir", dest='output_dir', action='store', default='data', help="Output directory")
parser.add_option('-e', "--ensembl", default=None, dest='ensembl', action='store', help="release version")
parser.add_option("--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")

(options, args) = parser.parse_args()

options.select = False

options.version = version

main.run(options)
