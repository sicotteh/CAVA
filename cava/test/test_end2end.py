import unittest
from cava.utils.data import Ensembl, Reference
from cava.utils import core
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
            raise Exception("Could not unzip genome: "+str(e)+"\n")

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
        #cls.genelist = ['BRCA2', 'BRCA1', "EPCAM", "MEN1", "SDHB", "PCSK9", "MUTYH"]
        #cls.transcriptlist = ['NM_000059.4','NM_007294.4', 'NM_002354.3', 'NM_001370259.2', 'NM_003000.3', 'NM_174936.4','NM_001048174.2']
        cls.codon_usage = ['1']
        cls.genelist=[]
        cls.transcriptlist = []
        cls.options = Options()
        cls.reference = Reference(cls.options)
        cls.ensembl = Ensembl(cls.options, cls.genelist, cls.transcriptlist, cls.codon_usage[0], cls.reference)


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
        # Non-Start Metionine mutation
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
        self.assertEqual('c.2389_2390del_p.Glu797ThrfsTer3', rec.variants[0].getFlag('CSN'))

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
        self.assertEqual('c.5558_5561dup_p.Ile1855ProfsTer26', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simple3utr(self):
        line = "17\t43044315\t438907\tT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*1363A>T', rec.variants[0].getFlag('CSN'))
    # Could be a dup or repeat
    # not within coding region, so not limited to %3
    #
    def test_hg38_negStrand_dup3utr(self):
        line = "17\t43044804\t803399\tC\tCTTTTTTTTT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*873A[19]%3B[28]', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_del3utr(self):
        line = "17\t43044804\t323411\tCTT\tC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.*873A[19]%3B[17]', rec.variants[0].getFlag('CSN'))

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
        self.assertEqual('c.5468-1G>T_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceAcceptor(self):
        line = "17\t43045749\t267605\tT\tTGTCCAACACCCACTCTCGGGTCACCACAGGTGCCTCACACATCTGCCCAATTGCTGGAGACAGA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-11_5520dup_p.Ala1843SerfsTer8', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tAA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5468-1_5468delinsTT_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceAcceptor(self):
        line = "17\t43045802\t632610\tGC\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
 #       self.assertEqual('c.5468-1del_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.5468-1del_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_simpleSpliceDonor(self):
        line = "17\t43047641\t125851\tA\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+2T>C_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
#        self.assertEqual('c.5467+1del_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.5467+1del_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_delinsSpliceDonor(self):
        line = "17\t43047641\t267604\tAC\tGG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5467+1_5467+2delinsCC_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_negStrand_dupSpliceDonor(self):
        line = "17\t43063871\t55424\tT\tTA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.5152+2dup_p.?', rec.variants[0].getFlag('CSN'))

    def test_hg38_allow_overlap_with_ends(self):
        line = "2\t47369285\toverlap_ends\tGGTGTGCGCTCCGCCCCGCCGCGCGCACAGAGCGCTAGTCCTTC\tGA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        csn=rec.variants[0].getFlag('CSN')
        self.assertEqual('5PU', rec.variants[0].getFlag('CLASS'))
        self.assertEqual(csn,"c.-221_-178delinsA","-221  is past the beginning of the transcript (-197) .. but that's OK")

    def test_hg38_snp_in_Met2(self):
        line = "13\t32316462\tsnp_in_met2\tT\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.2T>C_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_snp_in_Met1(self):
        line = "13\t32316461\tsnp_in_met1\tA\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.1A>C_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_snp_in_Met3(self):
        line = "13\t32316463\tsnp_in_met3\tG\tC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3G>C_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))

    def test_hg38_del_in_Met2(self):
        line = "13\t32316461\tindel_in_met2\tAT\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')  # So we can check value with debugger
        self.assertEqual('c.2del_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('1', rec.variants[0].getFlag('PROTPOS'))

    def test_hg38_del_in_Met1(self):
        # if  you shift the variants the right-most position, the  deleted A is the "A" in ATG of the CDS
        line = "13\t32316460\tindel_in_met1\tAA\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.1del_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('1', rec.variants[0].getFlag('PROTPOS'))
        self.assertEqual('frameshift_variant', rec.variants[0].getFlag('SO'))

    def test_hg38_ins_in_Met12(self):
        line = "13\t32316461\tins_in_met12\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	NM_000059.4(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.1_2insC_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('IM', rec.variants[0].getFlag('CLASS'))


    def test_hg38_ins_in_Met0_1(self):
        line = "13\t32316460\tins_in_met01\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        # NC_000013.11:g.32316462T>C	c(BRCA2):c.2T>C NP_000050.3:p.(?) p/M>T IM
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        rec.variants[0].getFlag('CSN')
        self.assertEqual('c.-1_1insC_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('5PU', rec.variants[0].getFlag('CLASS'))

    def test_hg38_ins_in_Ter(self): # Cannot be a Frameshift, needs to be a Stop Loss.
        line = "13\t32398769\tins_in_ter\tA\tAC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.10256_10257insC_p.Ter3419TyrextTer18', rec.variants[0].getFlag('CSN'))
        self.assertEqual('SL', rec.variants[0].getFlag('CLASS'))



    def test_hg38_splice_indel_donor(self): # 3' shifting at DNA level should work for exon/intron 3' shifting spf deletion at splice donor (BRCA2)
        line = "13\t32316527\tdonor_deletion\tGG\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.67+1del_p.?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('ESS', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('splice_donor_variant', rec.variants[0].getFlag('SO'))
        self.assertEqual('NC_000013.11:g.32316527G[2]%3B[1]', rec.variants[0].getFlag('HGVSg'))


    def test_hg38_splice_indel_exception_acceptor(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "13\t32325075\tacceptor_deletion\tGG\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.317del_p.Gly106GlufsTer15', rec.variants[0].getFlag('CSN'))
        self.assertEqual('FS', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('splice_region_variant|frameshift_variant', rec.variants[0].getFlag('SO'))
        self.assertEqual('NC_000013.11:g.32325075G[2]%3B[1]', rec.variants[0].getFlag('HGVSg'))

    def test_repeats(self):
        seq="TTTTT"
        [rep,nrep]=find_repeat_unit(seq)
        self.assertEqual("T",rep)
        self.assertEqual(5,nrep)

        seq = "TCTCT"
        [rep, nrep] = find_repeat_unit(seq)
        self.assertEqual("TCTCT", rep)
        self.assertEqual(1, nrep)


        seq = "TCTCTC"
        [rep, nrep] = find_repeat_unit(seq)
        self.assertEqual("TC", rep)
        self.assertEqual(3, nrep)

        seq = "GTCGTCGTC"
        [rep, nrep] = find_repeat_unit(seq)
        self.assertEqual("GTC", rep)
        self.assertEqual(3, nrep)

    def test_scan_for_repeat(self):
#        line = "13\t32329810\tinsertion_TA\tA\tATA\t.\t.\t.\tGT\t0/1\n"
        line = "13\t32329810\tdeletion_TA_\tATA\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        transcript = self.ensembl.findTranscripts(rec.variants[0])
        [dleftn,drightn,dfulln] = scan_for_repeat(rec.variants[0], self.reference)


        self.assertEqual(drightn[0], 32329817)  # deletion
        self.assertEqual(drightn[1], 32329807)  # deletion left shifted for HGVS .. position of first repeat
        self.assertEqual(drightn[2], 6)  # Number of ref
        self.assertEqual(drightn[3], 5)  # number of alt
        self.assertEqual(drightn[4], 'TA')
        self.assertEqual(drightn[5], 'TA')
        self.assertEqual(drightn[6], '')

        self.assertEqual(dleftn[0], 32329807)  # deletion
        self.assertEqual(dleftn[1], 32329807)  # deletion left shifted for HGVS .. position of first repeat
        self.assertEqual(dleftn[2], 6)  # Number of ref
        self.assertEqual(dleftn[3], 5)  # number of alt
        self.assertEqual(dleftn[4], 'TA')
        self.assertEqual(dleftn[5], 'TA')
        self.assertEqual(dleftn[6], '')

        self.assertEqual(dfulln[0], 32329807)  # deletion
        self.assertEqual(dfulln[1], 32329820)  # deletionleft shifted for HGVS .. position of first repeat
        self.assertEqual(dfulln[2], 'TATATATATATA')  # ref
        self.assertEqual(dfulln[3], 'TATATATATA')  # alt

        line = "13\t32329810\tdinsertions_TA_\tA\tATATA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        transcript = self.ensembl.findTranscripts(rec.variants[0])
        [ileftn, irightn, ifulln] = scan_for_repeat(rec.variants[0],  self.reference)
        self.assertEqual(irightn[0], 32329815)  # insertion, 1 bp after insertion site
        self.assertEqual(irightn[1], 32329807)  # insertion left shifted for HGVS .. position of first repeat
        self.assertEqual(irightn[2], 6)  # Number of ref
        self.assertEqual(irightn[3], 8)   # number of alt
        self.assertEqual(irightn[4], 'TA')
        self.assertEqual(irightn[5], '')
        self.assertEqual(irightn[6], 'TATA')

        self.assertEqual(ileftn[0], 32329807)  # insertion, 1 bp after insertion site
        self.assertEqual(ileftn[1], 32329807)  # insertion left shifted for HGVS .. position of first repeat
        self.assertEqual(ileftn[2], 6)  # Number of ref
        self.assertEqual(ileftn[3], 8)  # number of alt
        self.assertEqual(ileftn[4], 'TA')
        self.assertEqual(ileftn[5], '')
        self.assertEqual(ileftn[6], 'TATA')

        self.assertEqual(ifulln[0], 32329807)  # insertion, 1 bp after insertion site
        self.assertEqual(ifulln[1], 32329820)  # insertion left shifted for HGVS .. position of first repeat
        self.assertEqual(ifulln[2], 'TATATATATATA')  # ref
        self.assertEqual(ifulln[3], 'TATATATATATATATA')  # alt


    # BRCA1 is on the minus strand, so this variant should be left shifter and not be on the transcript.
    def test_hhg38_del_brca1_dupEnd(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "17\t43044291\tdel_TSS_shifted_into\tATGG\tA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NM_007294.4', rec.variants[0].getFlag('TRANSCRIPT'))
        self.assertEqual('Deletion', rec.variants[0].getFlag('TYPE'))

    # BRCA1 is on the minus strand, so this variant should be left shifter and not be on the transcript.
    def test_hhg38_ins_brca1_dupTSS(self): # 3' shifting at DNA level should shift Inside transcript
        line = "17\t43125364\tins_TSS_shifted_into\tC\tCGC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NM_007294.4', rec.variants[0].getFlag('TRANSCRIPT'))
        self.assertEqual('Insertion', rec.variants[0].getFlag('TYPE'))

    def test_hg38_snp(self):  # SNP in SELENON gene, after a UGA codon .. will error out if trim UGA codons (should be stop)
        line = "1\t25805147\tsbp_in_selenocysteine_code_gene\tA\tG\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NM_020451.3', rec.variants[0].getFlag('TRANSCRIPT'))
        self.assertEqual('Substitution', rec.variants[0].getFlag('TYPE'))


    def test_mito_ins_beg(self):  # Old code would fail, this close to the edge.
        line = "MT\t7\tmito_late_ins\tA\tACA\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual(4,rec.variants[0].alignonminus.pos)
        self.assertEqual('MT:3:T/TCA', rec.variants[0].alignonminus.id)
        # Normalization from the middle
        line = "MT\t6\tmito_mid_ins\tC\tCAC\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual(4, rec.variants[0].alignonminus.pos)
        self.assertEqual('MT:3:T/TCA', rec.variants[0].alignonminus.id)

#shifting at the ends did not use to work.
    def test_mito_del_end_with_padding_after(self):  # VCF format edge case
        line = "MT\t1\tmito_late_ins\tGAT\tT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual(1, rec.variants[0].alignonminus.pos)
        self.assertEqual('GA', rec.variants[0].alignonminus.ref)
        self.assertEqual('MT:1:GAT/T', rec.variants[0].alignonminus.id)

    def test_mito_ins_end(self):  # VCF format edge case
        line = "MT\t16568\tmito_late_ins\tT\tTGT\t.\t.\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual(16568, rec.variants[0].alignonminus.pos)
        self.assertEqual('TG', rec.variants[0].alignonminus.alt)
        self.assertEqual('MT:16567:A/ATG', rec.variants[0].alignonminus.id)
        self.assertEqual('', rec.variants[0].alignonplus.ref)
        self.assertEqual(16570, rec.variants[0].alignonplus.pos)  # Position points to 1 bp aafter insertion (longer than chromosome).. this will never be reported to user.
        self.assertEqual('MT:16569:G/GTG', rec.variants[0].alignonplus.id)


    #  # NM_015627.3(LDLRAP1): c.506_507delTT_p.Phe169Ter
    def test_hhg38_LDPRAP1_StopGainNotFS(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25562689\tStopgainFS\tTTT\tT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NM_015627.3', rec.variants[0].getFlag('TRANSCRIPT'))
        self.assertEqual('c.506_507del_p.Phe169Ter', rec.variants[0].getFlag('CSN'))
        self.assertEqual('SG', rec.variants[0].getFlag('CLASS'))
        self.assertEqual('frameshift_variant', rec.variants[0].getFlag('SO'))
        rec.output('VCF', None, self.options, self.genelist, self.transcriptlist, snplist=list(), stdout=True)

        # Add Test case for variant without a padded base (as occurs when splitting a multi-allele VCF

    def test_hhg38_LDPRAP1_not_dup(
            self):  # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563140\tnotdup_after_5C\tC\tCC\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
#        self.assertEqual('NM_015627.3(LDLRAP1):c.616dupT', rec.variants[0].getFlag('CSN'))

    def test_hhg38_LDPRAP1_not_dup2(
            self):  # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563140\tnotdup_after_5C\tC\tCCC\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
 #       self.assertEqual('NM_015627.3(LDLRAP1):c.616dupT', rec.variants[0].getFlag('CSN'))

    def test_hhg38_LDPRAP1_not_dup4(
            self):  # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563140\tnotdup_after_5C\tC\tCCCCC\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)

    def test_hhg38_LDPRAP1_repeat_del(
            self):  # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563138\tnotdup_after_5C\tCCC\tC\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)

    #       self.assertEqual('NM_015627.3(LDLRAP1):c.616dupT', rec.variants[0].getFlag('CSN'))

    # Add Test case for variant without a padded base (as occurs when splitting a multi-allele VCF
    def test_hhg38_LDPRAP1_dupT_last_before_splice(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563153\tdup_at_splice\tT\tTT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.616dup_p.?', rec.variants[0].getFlag('CSN'))


# insertion right at the last base of an exon .. at the junction.
# should be described as the dup of the last base of the exon (even if there is another 'T' on the 1st base of next exon
#chr1-25563153-T-TT
#OLD=NM_015627.3(LDLRAP1):c.616dupT
#NEW: NM_015627.3(LDLRAP1):c.616dupT


# ID=chr1-25563660-G-GAG
# OLD=NM_015627.3(LDLRAP1):c.617-2_617-1dupAG

    def test_hhg38_LDPRAP1_dupT_UTR(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "1\t25563660\tdup_acceptor\tG\tGAG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.617-2_617-1dup_p.?', rec.variants[0].getFlag('CSN'))

# ID=chr2-47805601-A-AT
# NM_000179.3(MSH6):c.3557-4dupT
    def test_hhg38_MSH6_dupT_intronic(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "2\t47805601\tdupT_intronic\tA\tAT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3557-16T[13]%3B[14]', rec.variants[0].getFlag('CSN'))


# chr2-47805710-A-ATA
# NM_000179.3(MSH6):c.3646+2_3646+3dupTA
# NM_000179.3(MSH6):c.3646+3_3646+4insTA
    def test_hhg38_MSH6_dupTA_intronic(self): # 3' shifting at DNA level should properly shift out a deletion of splice donor (BRCA2)
        line = "2\t47805710\tdupTA_intronic\tA\tATA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3646+2_3646+3dup_p.?', rec.variants[0].getFlag('CSN'))

# chr2-47806203-G-GAG
# OLD: NM_000179.3(MSH6):c.3647-2_3647-1dupAG
#   incorrect


# Already two 'G' there, so could be called dup, but because we are inserting two G ..it could be repeated allele
    # ... but the G's are inserted in the coding region, so it should NOT be a repeat (not a length 3 repeat)
    def test_hhg38_MSH6_dupGG_right_at_boundary_repeat(self):
        line = "2\t47806203\tdupG_junction\tG\tGGG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3647-1_3647dup_p.Arg1217GlufsTer12', rec.variants[0].getFlag('CSN'))
        self.assertEqual('NC_000002.12:g.47806203G[2]%3B[4]',rec.variants[0].getFlag('HGVSg'))
#        self.assertEqual('c.3647dup_p.Arg1217LysfsTer7', rec.variants[0].getFlag('CSN'))

# This GG cannot be described as a repeat, must be a dup
    def test_hhg38_MSH6_repeatTT_intron(self):
        line = "2\t47806192\tdupTT_intron\tT\tTTT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3647-14T[4]%3B[6]', rec.variants[0].getFlag('CSN'))

# This TTT should NOT be described as a repeat, even though it's in the coding region because 'T' is not a multiple of 3
# .. and although 3 T is a multiple of 3, other variants of 'T' repeat may not be.
#
    def test_hhg38_MSH6_repeatTTT_coding(self):
        line = "2\t47806222\tdupGG_junction\tT\tTTTT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3664_3666dup_p.Phe1222dup', rec.variants[0].getFlag('CSN'))

        # This GTT should be described as a repeat
    def test_hhg38_MSH6_repeatGTTGTT_coding(self):
        line = "2\t47804943\tdupGG_junction\tT\tTGTTGTT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.3470GTT[2]%3B[4]_p.Cys1158[1]%3B[3]', rec.variants[0].getFlag('CSN'))
        self.assertEqual('NC_000002.12:g.47804941GTT[2]%3B[4]',rec.variants[0].getFlag('HGVSg'))
# chr13-32319070-T-A,TA
# NM_000059.4(BRCA2):c.68-7T
# NM_000059.4(BRCA2):c.68-7delins A,NM_000059.4(BRCA2):c.68-6_68-5ins&quot;A

    def test_hhg38_delins(self):
        line = "13\t32319070\tdelins_\tT\tA,TA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.68-7T>A', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.68-6A[3]%3B[4]', rec.variants[1].getFlag('CSN')) # not in ESS or donor/acceptor, so no predicted impact on protein .. and variant not in protein.

# ID=chr13-32363178-G-GG
# NM_000059.4(BRCA2):c.7977-1dupG
# NM_000059.4(BRCA2):c.7977-1_7977insG
    def test_hhg38_dup(self):
        line = "13\t32363178\tdelins_before_splice_boundary\tG\tGG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.7977-1dup_p.?', rec.variants[0].getFlag('CSN'))



# chr19-11120523-G-GG
# NM_000527.5(LDLR):c.2140+1dupG
# sequence if G|GT  G|GGT .... rightmost sequetronnce is in in

    def test_hhg38_dup_donor_G_in_GT(self):
        line = "19\t11120523\tdup_donor_G_in_GT\tG\tGG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.2140+1dup_p.?', rec.variants[0].getFlag('CSN'))

#t
# Cannot be a repeat because after shifting to the repeat range,one end in tthe CDS .. and that requires repeat multiple of 3
#
    def test_hhg38_trip_donor_LDLR(self):
        line = "19\t11120523\tdup_donor_G_in_GT\tG\tGGG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
#        self.assertEqual('c.2140_2140+1G[2]%3B[4]', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.2140_2140+1dup_p.?', rec.variants[0].getFlag('CSN'))

#ID=chr2-47806858-T-TT   NM_000179.3(MSH6):c.4081dupT_p.Ter1361LeuextX3  c.4081dupT_p.Ter1361LeuextTer2  NM_000179.3(MSH6):c.4081dupT    NP_000170.1:p.(Ter1361LeuextTer2)   ---   >>NM_000179.3(MSH6):c.4081dupT_p.Ter1361LeuextTer2<<
    def test_hhg38_test_shortttest_ext(self):
        line = "2\t47806858\tshortttest_ext\tT\tTT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.4081dup_p.Ter1361LeuextTer2', rec.variants[0].getFlag('CSN'))

#ID=chr17-43124094-C-CCCT        NM_007294.4(BRCA1):c.2_3insAGG_p.?      c.2_3insAGG_p.Met1?     NM_007294.4(BRCA1):c.2_3insAGG  NP_009225.1:p.(Met1?)   ---   >>NM_007294.4(BRCA1):c.2_3insAGG_p.Met1?<<

    def test_hhg38_test_ins_inframe_firstcodon(self):
        line = "chr17\t43124094\tins_inframe_firstcodon\tC\tCCCT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.2_3insAGG_p.Met1?', rec.variants[0].getFlag('CSN'))

    def test_hhg38_inv(self):
        line = "chr17\t43124094\tinversion_minus_\tCATT\tAATG\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.-1_3inv_p.?', rec.variants[0].getFlag('CSN'))

# test insertion of sequence that is an inversion of sequence just 5' of insertion site (this example on minus strand)
    def test_hhg38_insinv(self):
        line = "chr17\t43124094\tinsinv\tC\tCAAT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.2_3ins-1_2inv_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('NC_000017.11:g.43124095_43124096insATA', rec.variants[0].getFlag('HGVSg'))  # Genomic is shifted 3' (right shifted by 1 position)


 # test insertion of sequence that is an inversion of sequence just 5' of insertion site (this example on minus strand)
 # it's only an inversion on the genomic, not on the cDNA.
    def test_hhg38_insinv_genomic(self):
        line = "chr17\t43124095\tinsinv_geno\tA\tATGG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.1_2insCCA_p.Met1?', rec.variants[0].getFlag('CSN'))
        self.assertEqual('NC_000017.11:g.43124096_43124097insGGT', rec.variants[0].getFlag('HGVSg')) # right shifted, so the 'inv' no longer applies

 # it's only an inversion on the genomic, not on the cDNA.
    def test_hhg38_notadup_bug(self):
        line = "chr17\t1340746\tnotadup\tG\tGCA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000017.11:g.1340746_1340747ins1340745_1340746inv', rec.variants[0].getFlag('HGVSg'))

    def test_trim_lead_sequence(self):
        line = "chr1\t25557155\tnotadup\tTCTCC\tTATCA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000001.11:g.25557156_25557159delinsATCA', rec.variants[0].getFlag('HGVSg'))

#    ID = chr3 - 37050680 - G - GGAGAGA
    def test_insert_repeat_no_repeat_on_ref(self):
        line = "chr3\t37050680\tnotarepeat\tG\tGGAGAGA\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('c.*27_*28insGAGAGA', rec.variants[0].getFlag('CSN'))


 # This deletion should be both repeat in
    def test_deletion_repeat_no_repeat_on_ref(self):
        line = "chr1\t25562689\tisrepeat\tTTT\tT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000001.11:g.25562689T[3]%3B[1]', rec.variants[0].getFlag('HGVSg'))
#        self.assertEqual('c.505_507T[3]%3B[1]_p.Phe169Ter', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.506_507del_p.Phe169Ter', rec.variants[0].getFlag('CSN')) # Not a repeat because repeat unit is not a multiple of 3 in coding region

    # Multiple transript overlap
#chr2	47806908	 A	ATTCAGACAACATTATGATCT
    def test_multiple_transcripts(self):
        line = "chr2\t47806908\tmulti_transcripts\tA\tATTCAGACAACATTATGATCT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000002.12:g.47806909_47806928dup', rec.variants[0].getFlag('HGVSg')) # should only be one value
        self.assertEqual('c.*49_*68dup:.',rec.variants[0].getFlag('CSN'))   # Two values
        rec.output('VCF', None, self.options, self.genelist, self.transcriptlist, snplist=list(), stdout=True) # Here we manuallyt he output checked

# chr2	47806751	35	CTT	C,CT
# Two alleles (Multiple Alleles by ','
#  so split by ';' then by ','. then by ':' (as long as no ) before )
# CAVA_CSN=c.4002-11_4002-10delTT,c.4002-10delT
# CAVA_HGVSc=.,.;CAVA_HGVSp=.,.;;CAVA_TYPE=Deletion,Deletion;
# CAVA_TRANSCRIPT=NM_000179.3,NM_000179.3;CAVA_GENE=MSH6,MSH6;CAVA_GENEID=MSH6,MSH6;CAVA_TRINFO=+/23810bp/10/4265bp/1360,+/23810bp/10/4265bp/1360;
# CAVA_LOC=In9/10,In9/10;CAVA_CSN=c.4002-11_4002-10del,c.4002-10del;CAVA_PROTPOS=.,.;CAVA_PROTREF=.,.;CAVA_PROTALT=.,.;
# CAVA_CLASS=INT,INT;CAVA_SO=intron_variant,intron_variant;CAVA_IMPACT=3,3;CAVA_ALTANN=c.4002-27_4002-26del,c.4002-27del;
# CAVA_ALTCLASS=.,.;CAVA_ALTSO=.,.;CAVA_ALTFLAG=AnnNotClassNotSO,AnnNotClassNotSO;CAVA_HGVSG=NC_000002.12:g.47806752T[18]%3B[16],NC_000002.12:g.47806769del;
# CAVA_HGVSc=NC_000002.12(NM_000179.3):c.4002-11_4002-10del,NC_000002.12(NM_000179.3):c.4002-10del;CAVA_HGVSp=.,.


    def test_multiple_alleles(self):
        line = "chr2\t47806751\tmulti_alleles\tCTT\tC,CT,<NON_REF>\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000002.12:g.47806752T[18]%3B[16]',
                         rec.variants[0].getFlag('HGVSg'))  # should only be one value
        # purely intronic variant, no protein and repeat-style annotation.
        self.assertEqual('c.4002-27T[18]%3B[16]', rec.variants[0].getFlag('CSN'))
        self.assertEqual('c.4002-27T[18]%3B[17]', rec.variants[1].getFlag('CSN'))
        # We manually checked the output using the debugger.. doesn't quite work as a unit test.
        rec.output('VCF', None, self.options, self.genelist, self.transcriptlist, snplist=list(), stdout=True)

#chr1	933741	.	T	TG	.
# IMP	CAVA_TYPE=Insertion;CAVA_TRANSCRIPT=NM_001385641.1;CAVA_GENE=SAMD11;CAVA_GENEID=SAMD11;CAVA_TRINFO=+/20653bp/14/3465bp/844;CAVA_LOC=In4/5;
# CAVA_CSN=c.843-2030_843-2023G[8]%3B[9];CAVA_PROTPOS=.;CAVA_PROTREF=.;CAVA_PROTALT=.;CAVA_CLASS=INT;CAVA_SO=intron_variant;CAVA_IMPACT=3;
# CAVA_ALTANN=c.843-2031_843-2030insG;CAVA_ALTCLASS=.;CAVA_ALTSO=.;CAVA_ALTFLAG=AnnNotClassNotSO;
# CAVA_HGVSg=NC_000001.11:g.933749dup;CAVA_HGVSc=NC_000001.11(NM_001385641.1):c.843-2030_843-2023G[8]%3B[9];CAVA_HGVSp=

    def test_repeat_in_intron(self):
        line = "chr1\t933741\trepeat_in_both\tT\tTG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000001.11:g.933742G[8]%3B[9]',
                         rec.variants[0].getFlag('HGVSg'))  # should only be one value
        # purely intronic variant, no protein and repeat-style annotation.
        self.assertEqual('c.843-2030G[8]%3B[9]', rec.variants[0].getFlag('CSN'))

#chr1	933548	.	A	AG	.	IMP	CAVA_TYPE=Insertion;CAVA_TRANSCRIPT=NM_001385641.1;CAVA_GENE=SAMD11;
# CAVA_GENEID=SAMD11;CAVA_TRINFO=+/20653bp/14/3465bp/844;CAVA_LOC=In4/5;
# CAVA_CSN=c.843-2223_843-2218G[6]%3B[7];CAVA_PROTPOS=.;CAVA_PROTREF=.;CAVA_PROTALT=.;
# CAVA_CLASS=INT;CAVA_SO=intron_variant;CAVA_IMPACT=3;CAVA_ALTANN=c.843-2224_843-2223insG;
# CAVA_ALTCLASS=.;CAVA_ALTSO=.;CAVA_ALTFLAG=AnnNotClassNotSO;
# CAVA_HGVSg=NC_000001.11:g.933554dup;CAVA_HGVSc=NC_000001.11(NM_001385641.1):c.843-2223_843-2218G[6]%3B[7];


    def test_repeat_in_intron(self):
        line = "chr1\t933548\trepeat_insertion_in_both\tT\tTG\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000001.11:g.933549G[6]%3B[7]',
                         rec.variants[0].getFlag('HGVSg'))  # should only be one value
        # purely intronic variant, no protein and repeat-style annotation.
        self.assertEqual('c.843-2223G[6]%3B[7]', rec.variants[0].getFlag('CSN'))


#chr1	1103993	.	C	CT	.	IMP	CAVA_TYPE=Insertion;CAVA_TRANSCRIPT=NM_017891.5;CAVA_GENE=C1orf159;CAVA_GENEID=C1orf159;
# CAVA_TRINFO=-/34268bp/10/1832bp/198;CAVA_LOC=In1/2;CAVA_CSN=c.-135-11891dup;CAVA_PROTPOS=.;CAVA_PROTREF=.;
# CAVA_PROTALT=.;CAVA_CLASS=5PU;CAVA_SO=5_prime_UTR_variant;CAVA_IMPACT=3;
# CAVA_ALTANN=c.-135-11901_-135-11900insA;CAVA_ALTCLASS=.;CAVA_ALTSO=.;
# CAVA_ALTFLAG=AnnNotClassNotSO;CAVA_HGVSg=NC_000001.11:g.1103994T[10]%3B[11];
# CAVA_HGVSc=;CAVA_HGVSp=.

    def test_repeat_in_intron(self):
        line = "chr1\t1103993\trepeat_insertion_in_both2\tC\tCT\t30\tPASS\t.\tGT\t0/1\n"
        rec = core.Record(line, self.options, None, self.reference)
        rec.annotate(self.ensembl, None, self.reference, None)
        self.assertEqual('NC_000001.11:g.1103994T[10]%3B[11]',
                         rec.variants[0].getFlag('HGVSg'))  # should only be one value
        # purely intronic variant, no protein and repeat-style annotation.
        self.assertEqual('c.-135-11891A[10]%3B[11]', rec.variants[0].getFlag('CSN'))


class Options:


    """Helper class for setting up testing options"""

    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        self.args = {#'ensembl': os.path.join(base_dir, 'data', 'RefSeq_small.gz'),
                     'ensembl': os.path.join(base_dir, 'data', 'cava_db.GRCh38.MANE.sorted.gz'),
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
