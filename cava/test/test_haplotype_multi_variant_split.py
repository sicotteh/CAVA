import os
import tempfile
import unittest
from collections import defaultdict
from urllib.parse import unquote

from cava.utils import core
from cava.utils import main


START_CODON_MULTI_CASES = [
    {
        'name': 'adjacent_met1_met2',
        'row': ('13', 32316461, '13_32316461_A_C;13_32316462_T_C', 'AT', 'CC'),
        'expected_split': ['c.1A>C_p.Met1?', 'c.2T>C_p.Met1?'],
    },
    {
        'name': 'slightly_separated_met1_met3',
        'row': ('13', 32316461, '13_32316461_A_C;13_32316463_G_C', 'ATG', 'CTC'),
        'expected_split': ['c.1A>C_p.Met1?', 'c.3G>C_p.Met1?'],
    },
    {
        'name': 'adjacent_met2_met3',
        'row': ('13', 32316462, '13_32316462_T_C;13_32316463_G_C', 'TG', 'CC'),
        'expected_split': ['c.2T>C_p.Met1?', 'c.3G>C_p.Met1?'],
    },
    {
        'name': 'three_adjacent_start_codon',
        'row': ('13', 32316461, '13_32316461_A_C;13_32316462_T_C;13_32316463_G_C', 'ATG', 'CCC'),
        'expected_split': ['c.1A>C_p.Met1?', 'c.2T>C_p.Met1?', 'c.3G>C_p.Met1?'],
    },
]


SERINE_SINGLE_EXPECTATIONS = {
    '17_43092111_A_C': 'c.3420T>G_p.Ser1140Arg',
    '17_43092111_A_G': 'c.3420T>C_p.Ser1140=',
    '17_43092111_A_T': 'c.3420T>A_p.Ser1140Arg',
    '17_43092112_C_A': 'c.3419G>T_p.Ser1140Ile',
    '17_43092112_C_G': 'c.3419G>C_p.Ser1140Thr',
    '17_43092112_C_T': 'c.3419G>A_p.Ser1140Asn',
    '17_43092113_T_A': 'c.3418A>T_p.Ser1140Cys',
    '17_43092113_T_C': 'c.3418A>G_p.Ser1140Gly',
    '17_43092113_T_G': 'c.3418A>C_p.Ser1140Arg',
    '17_43092114_A_C': 'c.3417T>G_p.Ser1139Arg',
    '17_43092114_A_G': 'c.3417T>C_p.Ser1139=',
    '17_43092114_A_T': 'c.3417T>A_p.Ser1139Arg',
    '17_43092115_C_A': 'c.3416G>T_p.Ser1139Ile',
    '17_43092115_C_G': 'c.3416G>C_p.Ser1139Thr',
    '17_43092115_C_T': 'c.3416G>A_p.Ser1139Asn',
    '17_43092116_T_A': 'c.3415A>T_p.Ser1139Cys',
    '17_43092116_T_C': 'c.3415A>G_p.Ser1139Gly',
    '17_43092116_T_G': 'c.3415A>C_p.Ser1139Arg',
}


class TestHaplotypeMultiVariantSplit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        cls.ref_path = os.path.join(cls.repo_root, 'cava', 'data', 'tmp.GRCh38.fa')
        cls.ens_path = os.path.join(cls.repo_root, 'cava', 'data', 'MANE.GRCh38.v1.1.refseq_genomic.db.gz')

    def _run_rows(self, rows):
        with tempfile.TemporaryDirectory() as td:
            cfg = os.path.join(td, 'cfg.txt')
            inp = os.path.join(td, 'in.vcf')
            outprefix = os.path.join(td, 'out')

            with open(cfg, 'w', encoding='utf-8') as f:
                f.write('@inputformat=VCF\n')
                f.write('@outputformat=VCF\n')
                f.write('@reference=' + self.ref_path + '\n')
                f.write('@ensembl=' + self.ens_path + '\n')
                f.write('@dbsnp=.\n')
                f.write('@logfile=FALSE\n')
                f.write('@prefix=FALSE\n')
                f.write('@chrom=.\n')

            with open(inp, 'w', encoding='utf-8') as f:
                f.write('##fileformat=VCFv4.2\n')
                f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                for chrom, pos, row_id, ref, alt in rows:
                    f.write(f'{chrom}\t{pos}\t{row_id}\t{ref}\t{alt}\t.\tPASS\t.\n')

            options = core.Options(cfg)
            options.args['parseHaplotype'] = True
            options.args['splitBasedOnProtein'] = True
            options.args['splitadjacentprotein'] = False

            with open(outprefix + '.vcf', 'w', encoding='utf-8') as out:
                core.writeHeader(options, '\n'.join(main.readHeader(inp)), out, False, '2.0.15')

            impactdir = {}
            for i, valuev in enumerate(options.args['impactdef'].split('|')):
                for c in valuev.split(','):
                    impactdir[c.strip()] = str(i + 1)

            class Copts:
                pass

            copts = Copts()
            copts.conf = cfg
            copts.input = inp
            copts.output = outprefix
            copts.stdout = False
            copts.threads = 1

            job = main.SingleJob(1, options, copts, 3, '', set(), set(), set(), impactdir, len(rows))
            job.run()

            outrows = []
            with open(outprefix + '.vcf', 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    pairs = {}
                    for item in cols[7].split(';'):
                        if '=' in item:
                            k, v = item.split('=', 1)
                            pairs[k] = v
                    outrows.append({
                        'chrom': cols[0],
                        'pos': int(cols[1]),
                        'id': cols[2],
                        'ref': cols[3],
                        'alt': cols[4],
                        'csn': unquote(pairs.get('CSN', '').split(':')[0]),
                        'orig': unquote(pairs.get('CAVA_ORIGHAPLOTYPE', '')),
                        'hap': unquote(pairs.get('CAVA_HAPLOTYPE', '')),
                    })
            return outrows

    def test_adjacent_and_slightly_separated_start_codon_haplotypes_split_to_singletons(self):
        # These atomic variants are already covered as single-variant tests in test_end2end.py.
        atoms = [
            ('13', 32316461, '13_32316461_A_C', 'A', 'C', 'c.1A>C_p.Met1?'),
            ('13', 32316462, '13_32316462_T_C', 'T', 'C', 'c.2T>C_p.Met1?'),
            ('13', 32316463, '13_32316463_G_C', 'G', 'C', 'c.3G>C_p.Met1?'),
        ]

        rows = [a[:5] for a in atoms] + [c['row'] for c in START_CODON_MULTI_CASES]
        outrows = self._run_rows(rows)

        by_orig = defaultdict(list)
        for r in outrows:
            if r['orig']:
                by_orig[r['orig']].append(r)

        for case in START_CODON_MULTI_CASES:
            records = by_orig[case['row'][2]]
            split_csns = [r['csn'] for r in records if r['hap'] != case['row'][2]]
            self.assertEqual(sorted(case['expected_split']), sorted(split_csns), msg=case['name'])

    def test_all_two_serine_snv_pairs_split_to_independent_expected_arrays(self):
        # Known consecutive Ser codons from existing single-variant region:
        # chr17:43092111-43092116 is ACTACT on genome (+), coding Ser1140/Ser1139 on transcript (-).
        chrom = '17'
        start = 43092111
        ref_seq = 'ACTACT'

        # Build all single SNVs across the 6 bp span.
        single_rows = []
        for i, ref_base in enumerate(ref_seq):
            pos = start + i
            for alt_base in 'ACGT':
                if alt_base == ref_base:
                    continue
                token = f'{chrom}_{pos}_{ref_base}_{alt_base}'
                single_rows.append((chrom, pos, token, ref_base, alt_base))

        # First pass: annotate singles to confirm the explicit expected array below matches
        # the already-existing single-variant tests for each atomic SNV.
        single_out = self._run_rows(single_rows)
        for r in single_out:
            if not r['orig']:
                self.assertEqual(SERINE_SINGLE_EXPECTATIONS[r['id']], r['csn'], msg=r['id'])

        # Group atoms by codon block (one SNV from each consecutive Ser codon).
        left_block = [r for r in single_rows if 43092111 <= r[1] <= 43092113]
        right_block = [r for r in single_rows if 43092114 <= r[1] <= 43092116]

        def build_haplotype_row(a, b):
            chrom_a, pos_a, token_a, _, alt_a = a
            chrom_b, pos_b, token_b, _, alt_b = b
            self.assertEqual(chrom_a, chrom_b)
            lo = min(pos_a, pos_b)
            hi = max(pos_a, pos_b)
            ref_span = list(ref_seq[lo - start: hi - start + 1])
            alt_span = ref_span[:]
            alt_span[pos_a - lo] = alt_a
            alt_span[pos_b - lo] = alt_b
            if pos_a < pos_b:
                row_id = token_a + ';' + token_b
                expected = [SERINE_SINGLE_EXPECTATIONS[token_a], SERINE_SINGLE_EXPECTATIONS[token_b]]
            else:
                row_id = token_b + ';' + token_a
                expected = [SERINE_SINGLE_EXPECTATIONS[token_b], SERINE_SINGLE_EXPECTATIONS[token_a]]
            return (chrom_a, lo, row_id, ''.join(ref_span), ''.join(alt_span)), expected

        multi_rows = []
        expected_arrays = []
        for a in left_block:
            for b in right_block:
                row, expected = build_haplotype_row(a, b)
                multi_rows.append(row)
                expected_arrays.append((row[2], expected))

        multi_out = self._run_rows(multi_rows)
        by_orig = defaultdict(list)
        for r in multi_out:
            by_orig[r['orig']].append(r)

        synonymous_seen = 0
        for row_id, expected in expected_arrays:
            records = by_orig[row_id]
            split_records = [r for r in records if r['hap'] != row_id]
            split_csns = [r['csn'] for r in split_records]
            self.assertEqual(sorted(expected), sorted(split_csns), msg=row_id)
            synonymous_seen += sum(1 for e in expected if '_p.Ser' in e and e.endswith('='))

        self.assertGreater(synonymous_seen, 0)


if __name__ == '__main__':
    unittest.main()
