#!/bin/env python3
import os
import sys
from optparse import OptionParser

from cava.utils  import main

with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

descr = 'CAVA (Clinical Annotation of VAriants) is a lightweight, fast and flexible NGS variant annotation tool that provides consistent transcript-level annotation. Limitation: does not call start gain variants in the UTR5'
epilog = '\nExample usage: python3 CAVA.py -c config.txt -i input.vcf -o output\n\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='python3 CAVA.py <options>'.format(version), version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--input", default='input.vcf', dest='input', action='store',
                  help="Input file name [default value: %default]")
parser.add_option('-o', "--output", default='output', dest='output', action='store',
                  help="Output file name prefix [default value: %default]")
parser.add_option('-c', "--config", default='config_template.txt', dest='conf', action='store',
                  help="Configuration file name [default value: %default]")
parser.add_option('-s', "--stdout", default=False, dest='stdout', action='store_true',
                  help="Write output to standard output [default value: %default]")
parser.add_option('-t', "--threads", default=1, dest='threads', action='store',
                  help="Number of threads [default value: %default]")
parser.add_option('--parseHaplotype', default=False, dest='parseHaplotype', action='store_true',
                  help='Parse semicolon-separated atomic haplotypes encoded in VCF ID [default value: %default]')
parser.add_option('--parseHaplotypee', default=False, dest='parseHaplotypee', action='store_true',
                  help='Deprecated alias for --parseHaplotype [default value: %default]')
parser.add_option('--splitBasedOnProtein', default=False, dest='splitBasedOnProtein', action='store_true',
                  help='Split parsed haplotypes into subsets and reannotate [default value: %default]')
parser.add_option('--splitadjacentprotein', default=False, dest='splitadjacentprotein', action='store_true',
                  help='Emit optional adjacent protein split outputs for parsed haplotypes [default value: %default]')

# Backward-compatible one-dash alias requested by haplotype spec.
argv = [('--splitBasedOnProtein' if x == '-splitBasedOnProtein' else x) for x in sys.argv]
(copts, args) = parser.parse_args(argv[1:])

if copts.parseHaplotypee:
    copts.parseHaplotype = True

main.run(copts, version)
