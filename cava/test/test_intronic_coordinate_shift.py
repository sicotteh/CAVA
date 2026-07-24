import unittest

from cava.utils.core import Transcript, Variant
from cava.utils import csn
from cava.utils.csn import calculateCSNCoordinates, transformToCSNCoordinate


class MockFasta:
    def __init__(self, seq):
        self.seq = seq
        self.references = ['1']
        self.lengths = [len(seq)]

    def get_reference_length(self, chrom):
        return len(self.seq)


class MockReference:
    def __init__(self, seq):
        self.fastafile = MockFasta(seq)
        self.reflens = {'1': len(seq), 'chr1': len(seq)}

    def getReference(self, chrom, start, end):
        if chrom not in ('1', 'chr1'):
            raise Exception('Unknown chromosome')
        if start < 1 or end > len(self.fastafile.seq):
            raise Exception('Out of range')
        return self.fastafile.seq[start - 1:end]


class TestIntronicCoordinateShift(unittest.TestCase):
    def setUp(self):
        # Exons are 100..110 and 121..130 (1-based inclusive), intron is 111..120.
        self.transcript = Transcript(
            '\t'.join(
                ['TX1', 'GENE', 'G1', '.', '1', '1', '99', '130', '1', '100', '130', '99', '110', '120', '130']
            )
        )
        # Intron is a run of A's so right-normalization can move indels across midpoint.
        seq = 'C' * 99 + 'G' * 11 + 'A' * 10 + 'T' * 10 + 'C' * 100
        self.reference = MockReference(seq)

    def test_deletion_switches_exon_anchor_after_right_shift(self):
        # VCF AA>A at 112 normalizes to a single-base deletion at 113.
        var = Variant('1', 112, 'AA', 'A')
        shifted = var.alignOnPlusStrand(self.reference)

        before = calculateCSNCoordinates(var, self.transcript)
        after = calculateCSNCoordinates(shifted, self.transcript)

        self.assertEqual(('11', 3), (before[0], before[1]))
        # After right-normalization, deletion moves to intron end and flips to downstream-exon anchor.
        self.assertEqual(('12', -1), (after[0], after[1]))

    def test_insertion_switches_exon_anchor_after_right_shift(self):
        var = Variant('1', 112, 'A', 'AA')
        shifted = var.alignOnPlusStrand(self.reference)

        before = calculateCSNCoordinates(var, self.transcript)
        after = calculateCSNCoordinates(shifted, self.transcript)

        self.assertEqual(('11', 2), (before[0], before[1]))
        self.assertEqual(('12', -1), (after[0], after[1]))

    def test_odd_intron_central_base_uses_plus_notation(self):
        # HGVS numbering rule: in odd-length introns, the central nucleotide is + relative to upstream exon.
        odd_intron_tx = Transcript(
            '\t'.join(
                ['TX2', 'GENE', 'G1', '.', '1', '1', '99', '130', '1', '100', '130', '99', '110', '119', '130']
            )
        )
        central_pos = 115  # Intron 111..119
        self.assertEqual(('11', 5, 0), transformToCSNCoordinate(central_pos, odd_intron_tx))

    def test_single_repeat_unit_deletion_uses_del_not_repeat(self):
        var = Variant('1', 112, 'AA', 'A')
        shifted = var.alignOnPlusStrand(self.reference)

        annotation, _, alt_annotation = csn.getAnnotation(shifted, self.transcript, self.reference, '', None)

        self.assertEqual('c.12-1del', annotation.getAsString())
        self.assertIsNone(alt_annotation)

    def test_single_repeat_unit_multibase_deletion_uses_del_not_repeat(self):
        # Exons 100..110 and 123..132, intron 111..122 contains six TA copies.
        transcript = Transcript(
            '\t'.join(
                ['TX3', 'GENE', 'G1', '.', '1', '1', '99', '132', '1', '100', '132', '99', '110', '122', '132']
            )
        )
        reference = MockReference('C' * 99 + 'G' * 11 + 'TA' * 6 + 'T' * 10 + 'C' * 100)
        var = Variant('1', 111, 'TAT', 'T')
        shifted = var.alignOnPlusStrand(reference)

        annotation, _, alt_annotation = csn.getAnnotation(shifted, transcript, reference, '', None)

        self.assertEqual('c.12-1_12del', annotation.getAsString())
        self.assertIsNone(alt_annotation)

    def test_single_repeat_unit_multibase_insertion_uses_dup(self):
        transcript = Transcript(
            '\t'.join(
                ['TX3', 'GENE', 'G1', '.', '1', '1', '99', '132', '1', '100', '132', '99', '110', '122', '132']
            )
        )
        reference = MockReference('C' * 99 + 'G' * 11 + 'TA' * 6 + 'T' * 10 + 'C' * 100)
        var = Variant('1', 112, 'T', 'TTA')
        shifted = var.alignOnPlusStrand(reference)

        annotation, _, alt_annotation = csn.getAnnotation(shifted, transcript, reference, '', None)

        self.assertEqual('c.12-1_12dup', annotation.getAsString())
        self.assertIsNone(alt_annotation)


if __name__ == '__main__':
    unittest.main()
