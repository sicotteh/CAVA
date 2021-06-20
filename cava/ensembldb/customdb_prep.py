import argparse
import sys
import textwrap
from os import path

sys.path.insert(0, path.dirname(path.realpath(__file__)) + '/pysamdir')
import pysam


# Use Tabix to index the custom database file
def indexFile(input_file):
    sys.stdout.write('Compressing file... ')
    sys.stdout.flush()
    pysam.tabix_compress(input_file, input_file + '.gz', force=True)
    sys.stdout.write('OK\n')
    sys.stdout.write('Indexing output file... ')
    sys.stdout.flush()
    pysam.tabix_index(input_file + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)
    sys.stdout.write('OK\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
     This utility is designed to compress a custom cava database from a pre-sorted input text file.
     
     The input file must be sorted on chromosome, and then by the transcript start position. 

     Internal file format is expected to be tab separated:
     col 1  -- Transcript title
     col 2  -- Gene name
     col 3  -- Gene id
     col 4  -- TRINFO
     col 5  -- Chromosome
     col 6  -- Strand
     col 7  -- Transcript start position
     col 8  -- Transcript end position
     col 9  -- Coding relative start position
     col 10 -- Coding start position
     col 11 -- Coding end position
     col 12 -- Exon 1 start position
     col 13 -- Exon 1 end position
     col 14 -- Exon 2 start position
     col 15 -- Exon 2 end position
     exons continue
         '''))
    parser.add_argument('-i', dest='input_file', required='true', help='path to the input file')
    args = parser.parse_args()

    if not path.isfile(args.input_file):
        print
        "This does not appear to be a file"

    indexFile(path.abspath(args.input_file))
