import unittest
from cava.utils.data import Ensembl, Reference
from cava.utils import core
import os
import pycurl
import sys


def check_materials():
    base_dir = os.path.dirname(os.path.dirname(__file__))

    if not os.path.exists(os.path.join(base_dir, 'data')):
        print('Making data directory: {}'.format(os.path.join(base_dir, 'data')))
        os.mkdir(os.path.join(base_dir, 'data'))

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.GRCh38.fa.fai')) and not os.path.exists(
            os.path.join(base_dir, 'data', 'tmp.GRCh38.fa')):
        print('Downloading build 38')
        c = pycurl.Curl()
        c.setopt(c.URL, 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')
        with open(os.path.join(base_dir, 'data', 'tmp.GRCh38.fa.gz'), 'wb') as f:
            c.setopt(c.WRITEDATA, f)
            c.perform()

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.GRCh38.fa')):
        try:
            print('Unzipping genome')
            os.system("bgzip -d " + os.path.join(base_dir, 'data', 'tmp.GRCh38.fa.gz'))
        except Exception as e:
            raise Exception("Could not unzip genome.")

    if not os.path.exists(os.path.join(base_dir, 'data', 'tmp.GRCh38.fa.fai')):
        try:
            print('Indexing genome')
            os.system("samtools faidx " + os.path.join(base_dir, 'data', 'tmp.GRCh38.fa'))
        except Exception as e:
            raise Exception("Unable to index the reference genome. Do you have samtools in your path?")


class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        check_materials()
        cls.genelist = ['BRCA2', 'BRCA1', "EPCAM", "MEN1", "SDHB", "PCSK9", "MUTYH"]
        cls.transcriptlist = ['NM_000059.4','NM_007294.4', 'NM_002354.3', 'NM_001370259.2', 'NM_003000.3', 'NM_174936.4','NM_001048174.2']
        cls.codon_usage = ['1']
        cls.options = Options()
        cls.ensembl = Ensembl(cls.options, cls.genelist, cls.transcriptlist, cls.codon_usage[0])
        cls.reference = Reference(cls.options)

    def test_hg38_negStrand_simple_nc(self):
        line = "17\t43045705\t869060\tT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5565A>T_p.Ile1855=', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simplemissense(self):
        line = "17\t43051071\tVCV000017694.14\tA\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5324T>G_p.Met1775Arg', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simplestop(self):
        line = "17\t43045711\t55630\tG\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5559C>G_p.Tyr1853Ter', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_fsdel(self):
        line = "17\t43093140\trs80357695\tTTC\tT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.2389_2390delGA_p.Glu797ThrfsTer3', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_fsdelins(self):
        line = "17\t43045707\trs483353103\tT\tGGATCC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5563delinsGGATCC_p.Ile1855GlyfsTer69', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_fsdup(self):
        line = "17\t43045708\t949924\tC\tCAGGT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5558_5561dupACCT_p.Ile1855ProfsTer26', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simple3utr(self):
        line = "17\t43044315\t438907\tT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*1363A>T', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dup3utr(self):
        line = "17\t43044804\t803399\tC\tCTTTTTTTTT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*865_*873dup9', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_del3utr(self):
        line = "17\t43044804\t323411\tCTT\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*872_*873delAA', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delins3utr(self):
        line = "17\t43044804\t323411\tCTT\tCG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*872_*873delinsC', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simpleSpliceAcceptor(self):
        line = "17\t43045803\t864908\tC\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1G>T', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceAcceptor(self):
        line = "17\t43045749\t267605\tT\tTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATTGCTGGAGACAGA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-11_5520dup64_p.Ala1843SerfsTer8', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tAA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1_5468delinsTT', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1delG', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simpleSpliceDonor(self):
        line = "17\t43047641\t125851\tA\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+2T>C', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+1delG', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tGG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+1_5467+2delinsCC', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceDonor(self):
        line = "17\t43063871\t55424\tT\tTA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5152+2dupT', rec.variants[0].getFlag('CSN'))

    def test_hg38_allow_overlap_with_ends(self):
        line = "2\t47369285\toverlap_ends\tGGTGTGCGCTCCGCCCCGCCGCGCGCACAGAGCGCTAGTCCTTC\tGA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        csn=rec.variants[0].getFlag('CSN')
        sys.stderr.write("CSN for large del overlapping ends="+csn+"\n")
        self.assertEqual('5PU', rec.variants[0].getFlag('CLASS'))
        self.assertEqual(csn,"c.-221_-178delinsA","-221  is past the beginning of the transcript (-197) .. but that's OK")

    def test_hg38_snp_in_Met2(self):
        line = "13\t32316462\tsnp_in_met2\tT\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.2T>C_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_snp_in_Met1(self):
        line = "13\t32316461\tsnp_in_met1\tA\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.1A>C_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_snp_in_Met3(self):
        line = "13\t32316463\tsnp_in_met3\tG\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3G>C_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_del_in_Met2(self):
        line = "13\t32316461\tindel_in_met2\tAT\tA\t.\t.\t.\tGT\t0/1\n"

        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')  # So we can check value with debugger
        self.assertEqual('c.2delT_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_del_in_Met1(self):
        line = "13\t32316460\tindel_in_met1\tAA\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.1delA_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('frameshift_variant', rec.variants[0].getFlag('SO'))

    def test_hg38_ins_in_Met12(self):
        line = "13\t32316461\tins_in_met12\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.1_2insC_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))


    def test_hg38_ins_in_Met0_1(self):
        line = "13\t32316460\tins_in_met01\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	c(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.-1_1insC', rec.variants[0].getFlag('CSN'))
        self.assertEqual('5PU', rec.variants[0].getFlag('CLASS'))

    def test_hg38_ins_in_Ter(self): # Cannot be a Frameshift, needs to be a Stop Loss.
        line = "13\t32398769\tins_in_ter\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	c(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.10256_10257insC_p.Ter3419TyrextX19', rec.variants[0].getFlag('CSN'))
        self.assertEqual('SL', rec.variants[0].getFlag('CLASS'))



    def test_hg38_splice_indel_donor(self): # 3' shifting at DNA level should work for exon/intron 3' shifting spf deletion at splice donor (BRCA2)
        line = "13\t32316526\tdonor_deletion\tGG\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.67+1delG', rec.variants[0].getFlag('CSN'))
        self.assertEqual('ESS', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('splice_donor_variant', rec.variants[0].getFlag('SO'))


    def test_hg38_splice_indel_exception_acceptor(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "13\t32325074\tacceptor_deletion\tGG\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.317delG_p.Gly106GlufsTer15', rec.variants[0].getFlag('CSN'))
        self.assertEqual('FS', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('splice_region_variant|frameshift_variant', rec.variants[0].getFlag('SO'))



class Options:
    """Helper class for setting up testing options"""

    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        self.args = {'ensembl': os.path.join(base_dir, 'data', 'RefSeq_small.gz'),
                     'logfile': None,
                     'reference': os.path.join(base_dir, 'data', 'tmp.GRCh38.fa'),
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
