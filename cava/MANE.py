from optparse import OptionParser
from cava.ensembldb import mane_db_prep as main
import os


def _split_csv_options(values):
    result = set()
    for value in values or []:
        for item in value.split(","):
            item = item.strip()
            if item:
                result.add(item)
    return result

with open(os.path.join(os.path.dirname(__file__), "VERSION")) as version_file:
    version = version_file.read().strip()

# Command line argument parsing
descr = "MANE.py" + version
epilog = (
    "\nExample usage: MANE.py -e 1.3 -D mane_1.3\n"
    "Note: by default, hg19 will be created using crossmap\n"
    "Version: {}\n"
    "\n".format(version)
)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(
    usage="\n\nMANE.py <options>", version=version, description=descr, epilog=epilog
)
parser.add_option(
    "-D",
    "--outdir",
    dest="output_dir",
    action="store",
    default="data",
    help="Output directory",
)
parser.add_option(
    "-e",
    "--mane_version",
    default=None,
    dest="ensembl",
    action="store",
    help="MANE release version",
)
parser.add_option(
    "--no_hg19",
    action="store_false",
    default=True,
    dest="no_hg19",
    help="Set this to skip hg19 builds",
)
parser.add_option(
    "--include-alt-gene",
    action="append",
    default=[],
    dest="include_alt_genes",
    help=(
        "Gene symbol to retain when MANE annotation is on an alternate contig. "
        "Repeat option or use commas. Example: --include-alt-gene GSTT1"
    ),
)
parser.add_option(
    "--include-alt-transcript",
    action="append",
    default=[],
    dest="include_alt_transcripts",
    help=(
        "Transcript accession to retain on a non-primary contig. Repeat option "
        "or use commas. A versioned accession also matches MANE suffix forms "
        "such as NM_000853.4_1."
    ),
)

options, args = parser.parse_args()

options.include_alt_genes = _split_csv_options(options.include_alt_genes)
options.include_alt_transcripts = _split_csv_options(options.include_alt_transcripts)

options.select = False

options.version = version
if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)

main.run(options)
