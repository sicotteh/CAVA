import unittest
from cava.utils.data import Ensembl, Reference
import cava.utils.core as core
from cava.utils.csn import find_repeat_unit
from cava.utils.csn import scan_for_repeat
import os
import pycurl
import sys

def check_materials():
    base_dir = os.path.dirname(os.path.dirname(__file__)) # This file in cava/test .. this points to base_dir=cava


    if not os.path.exists(os.path.join(base_dir, 'data')):
        print('Making data directory: {}'.format(os.path.join(base_dir, 'data')))
        os.mkdir(os.path.join(base_dir, 'data'))

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.hg19.fa')):
        print('Downloading build 19 fasta')
        c = pycurl.Curl()
        c.setopt(c.URL, 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz')
        with open(os.path.join(base_dir, 'data', 'tmp.hg19.fa.gz'), 'wb') as f:
            c.setopt(c.WRITEDATA, f)
            c.perform()

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.hg19.fa')):
        try:
            print('Unzipping genome')
            os.system("bgzip -d " + os.path.join(base_dir, 'data', 'tmp.hg19.fa.gz'))
        except Exception as e:
            raise Exception("Could not unzip genome: "+str(e)+"\n")

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.hg19.fa.fai')):
        try:
            print('Indexing genome')
            os.system("samtools faidx " + os.path.join(base_dir, 'data', 'tmp.hg19.fa'))
        except Exception as e:
            raise Exception("Unable to index the reference genome. Do you have samtools in your path?")


class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        check_materials()

        cls.codon_usage = ['1']
        cls.genelist=[]
        cls.transcriptlist = []
        cls.options = Options()
        cls.reference = Reference(cls.options)
        cls.ensembl = Ensembl(cls.options, cls.genelist, cls.transcriptlist, cls.codon_usage[0], cls.reference)

# 1:26131654:G/A
    def testProteinTooShort(self):  # GAA(NM_001079804.3):c.-32-13T>G
        line = "chr1\t26131654\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        a = 1 # anchor point for debugging.
        self.assertEqual('ESS', rec.variants[0].getFlag('CLASS'))

        self.assertEqual('NC_000013.11:g.32339700A[7]%3B[5]', rec.variants[0].getFlag('HGVSg'))

        self.assertEqual('c.1_3dup_p.Met1dup', rec.variants[0].getFlag('CSN'))

class Options:


    """Helper class for setting up testing options"""

    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        self.args = {
                     'ensembl': os.path.join(base_dir, 'data', 'ENST75_GRCh37.gz'),
                     'logfile': None,
                     'reference': os.path.join(base_dir, 'data', 'tmp.hg19.fa'),
                     'inputformat': 'VCF',
                     'type': 'ALL',
                     'ontology': 'BOTH',
                     'givealt': True,
                     'ssrange': 4,
                     'impactdef': 'SG,ESS,FS|SS5,IM,SL,EE,IF,NSY|SY,SS,INT,5PU,3PU',
                        'selenofile' : 'SECIS_in_refseq_pos.homo_sapiens.ensembl_GRCh37_75.txt'
                     }
        self.dbsnp = None
        self.transcript2protein = None
        super().__init__()
        self.transcript2protein = core.read_dict(self, 'transcript2protein')


if __name__ == '__main__':
    unittest.main()
