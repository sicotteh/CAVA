import unittest
from cava_.data import Ensembl, Reference
from cava_ import core

class MyTestCase(unittest.TestCase):

    def __init__(self):
        self.genelist = ['BRCA1']
        self.transcriptlist = ['NM_007294.4']
        self.codon_usage = ['1']
        self.options = Options()
        self.ensembl = Ensembl(options, genelist, transcriptlist, codon_usage[0])
        self.reference = Reference(options)

    def test_hg38_negStrand_simple_nc(self):
        line = "17\t43045705\t869060\tT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5565A>T_p.Ile1855=', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simplemissense(self):
        line = "17\t43051071\tVCV000017694.14\tA\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5324T>G_p.Met1775Arg', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simplestop(self):
        line = "17\t43045711\t55630\tG\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5559C>G_p.Tyr1853Ter', rec.variants[0].getFlag('CSN'))


    def test_hg38_negStrand_fsdel(self):
        line = "17\t43093140\trs80357695\tTTC\tT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.2389_2390delGA_p.Glu797ThrfsTer3', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_fsdelins(self):
        line = "17\t43045707\trs483353103\tT\tGGATCC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5563delinsGGATCC_p.Ile1855GlyfsTer69', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_fsdup(self):
        line = "17\t43045708\t949924\tC\tCAGGT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5558_5561dupACCT_p.Ile1855ProfsTer26', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simple3utr(self):
        line = "17\t43044315\t438907\tT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.+1363A>T', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dup3utr(self):
        line = "17\t43044804\t803399\tC\tCTTTTTTTTT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.+865_+873dup9', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_del3utr(self):
        line = "17\t43044804\t323411\tCTT\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.+872_+873delAA', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delins3utr(self):
        line = "17\t43044804\t323411\tCTT\tCG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.+872_+873delinsC', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simpleSpliceAcceptor(self):
        line = "17\t43045803\t864908\tC\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1G>T', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceAcceptor(self):
        line = "17\t43045749\t267605\tT\tTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATTGCTGGAGACAGA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-11_5520dup64_p.Ala1843SerfsTer8', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tAA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1_5468delinsTT', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1delG', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simpleSpliceDonor(self):
        line = "17\t43047641\t125851\tA\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+2T>C', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+1delG', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tGG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+1_5467+2delinsCC', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceDonor(self):
        line = "17\t43063871\t55424\tT\tTA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, options, None)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5152+2dupT', rec.variants[0].getFlag('CSN'))


class Options:
    """Helper class for setting up testing options"""
    def __init__(self):
        self.args = {'ensembl': 'data/RefSeq_small.gz',
                 'logfile': None,
                 'reference': 'data/tmp.GRCh38.fa',
                 'inputformat': 'VCF',
                 'type': 'ALL',
                 'ontology': 'BOTH',
                 'givealt': None,
                 'ssrange': 4,
                 'impactdef': 'SG,ESS,FS|SS5,IM,SL,EE,IF,NSY|SY,SS,INT,5PU,3PU'
                 }
        self.dbsnp = None
        self.transcript2protein = None
        super().__init__()
        self.transcript2protein = core.read_dict(self, 'transcript2protein')


if __name__ == '__main__':
    unittest.main()
