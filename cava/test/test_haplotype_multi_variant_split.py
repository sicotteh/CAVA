import io
import os
import tempfile
import unittest
from collections import defaultdict
from contextlib import redirect_stdout
from itertools import product
from unittest.mock import patch
from urllib.parse import unquote

from cava.utils import core
from cava.utils import haplotype
from cava.utils import main
from cava.utils.data import Ensembl, Reference

START_CODON_MULTI_CASES = [
    {
        "name": "adjacent_met1_met2",
        "row": ("13", 32316461, "13_32316461_A_C;13_32316462_T_C", "AT", "CC"),
        "expected_split": ["c.1A>C_p.Met1?", "c.2T>C_p.Met1?"],
    },
    {
        "name": "slightly_separated_met1_met3",
        "row": ("13", 32316461, "13_32316461_A_C;13_32316463_G_C", "ATG", "CTC"),
        "expected_split": ["c.1A>C_p.Met1?", "c.3G>C_p.Met1?"],
    },
    {
        "name": "adjacent_met2_met3",
        "row": ("13", 32316462, "13_32316462_T_C;13_32316463_G_C", "TG", "CC"),
        "expected_split": ["c.2T>C_p.Met1?", "c.3G>C_p.Met1?"],
    },
    {
        "name": "three_adjacent_start_codon",
        "row": (
            "13",
            32316461,
            "13_32316461_A_C;13_32316462_T_C;13_32316463_G_C",
            "ATG",
            "CCC",
        ),
        "expected_split": ["c.1A>C_p.Met1?", "c.2T>C_p.Met1?", "c.3G>C_p.Met1?"],
    },
]


SERINE_SINGLE_EXPECTATIONS = {
    "17_43092111_A_C": "c.3420T>G_p.Ser1140Arg",
    "17_43092111_A_G": "c.3420T>C_p.Ser1140=",
    "17_43092111_A_T": "c.3420T>A_p.Ser1140Arg",
    "17_43092112_C_A": "c.3419G>T_p.Ser1140Ile",
    "17_43092112_C_G": "c.3419G>C_p.Ser1140Thr",
    "17_43092112_C_T": "c.3419G>A_p.Ser1140Asn",
    "17_43092113_T_A": "c.3418A>T_p.Ser1140Cys",
    "17_43092113_T_C": "c.3418A>G_p.Ser1140Gly",
    "17_43092113_T_G": "c.3418A>C_p.Ser1140Arg",
    "17_43092114_A_C": "c.3417T>G_p.Ser1139Arg",
    "17_43092114_A_G": "c.3417T>C_p.Ser1139=",
    "17_43092114_A_T": "c.3417T>A_p.Ser1139Arg",
    "17_43092115_C_A": "c.3416G>T_p.Ser1139Ile",
    "17_43092115_C_G": "c.3416G>C_p.Ser1139Thr",
    "17_43092115_C_T": "c.3416G>A_p.Ser1139Asn",
    "17_43092116_T_A": "c.3415A>T_p.Ser1139Cys",
    "17_43092116_T_C": "c.3415A>G_p.Ser1139Gly",
    "17_43092116_T_G": "c.3415A>C_p.Ser1139Arg",
}


class TestHaplotypeMultiVariantSplit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        cls.ref_path = os.path.join(cls.repo_root, "cava", "data", "tmp.GRCh38.fa")
        cls.ens_path = os.path.join(
            cls.repo_root, "cava", "data", "MANE.GRCh38.v1.1.refseq_genomic.db.gz"
        )
        cls.reference = Reference(type("Opt", (), {"args": {"reference": cls.ref_path}})())

    @staticmethod
    def _splice_severity(class_value, so_value):
        classes = [x.strip() for x in str(class_value).split(":") if x.strip()]
        sos = [x.strip() for x in str(so_value).split(":") if x.strip()]

        for c in classes:
            if c == "ESS":
                return 2
        for s in sos:
            if (
                "splice_acceptor_variant" in s
                or "splice_donor_variant" in s
                or "splice_donor_5th_base_variant" in s
            ):
                return 2

        for c in classes:
            if c in {"SS", "SS5", "EE"}:
                return 1
        for s in sos:
            if "splice_region_variant" in s:
                return 1

        return 0

    def _build_row_from_tokens(self, tokens):
        atoms = [haplotype.parse_atomic_token(tok) for tok in tokens]
        pos, ref, alt, row_id = haplotype.build_subset_vcf_fields(
            self.reference, atoms[0].chrom, atoms
        )
        return (atoms[0].chrom, pos, row_id, ref, alt)

    def _run_rows(
        self, rows, split_based=True, split_adj=False, return_raw=False,
        load_all=True, output_format="VCF",
    ):
        with tempfile.TemporaryDirectory() as td:
            cfg = os.path.join(td, "cfg.txt")
            inp = os.path.join(td, "in.vcf")
            outprefix = os.path.join(td, "out")

            with open(cfg, "w", encoding="utf-8") as f:
                f.write("@inputformat=VCF\n")
                f.write("@outputformat=" + output_format + "\n")
                f.write("@reference=" + self.ref_path + "\n")
                f.write("@ensembl=" + self.ens_path + "\n")
                f.write("@dbsnp=.\n")
                f.write("@logfile=FALSE\n")
                f.write("@prefix=FALSE\n")
                f.write("@chrom=.\n")

            with open(inp, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for chrom, pos, row_id, ref, alt in rows:
                    f.write(f"{chrom}\t{pos}\t{row_id}\t{ref}\t{alt}\t.\tPASS\t.\n")

            options = core.Options(cfg)
            options.args["parseHaplotype"] = True
            options.args["splitBasedOnProtein"] = bool(split_based)
            options.args["splitadjacentprotein"] = bool(split_adj)
            options.args["loadalltranscripts"] = load_all

            output_path = outprefix + (".vcf" if output_format == "VCF" else ".txt")
            with open(output_path, "w", encoding="utf-8") as out:
                core.writeHeader(
                    options, "\n".join(main.readHeader(inp)), out, False, "2.0.15"
                )

            impactdir = {}
            for i, valuev in enumerate(options.args["impactdef"].split("|")):
                for c in valuev.split(","):
                    impactdir[c.strip()] = str(i + 1)

            class Copts:
                pass

            copts = Copts()
            copts.conf = cfg
            copts.input = inp
            copts.output = outprefix
            copts.stdout = False
            copts.threads = 1

            job = main.SingleJob(
                1, options, copts, 3, "", set(), set(), set(), impactdir, len(rows)
            )
            job.run()

            outrows = []
            with open(output_path, "r", encoding="utf-8") as f:
                if return_raw:
                    return f.read()
                for line in f:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    pairs = {}
                    for item in cols[7].split(";"):
                        if "=" in item:
                            k, v = item.split("=", 1)
                            pairs[k] = v
                    outrows.append(
                        {
                            "chrom": cols[0],
                            "pos": int(cols[1]),
                            "id": cols[2],
                            "ref": cols[3],
                            "alt": cols[4],
                            "csn": unquote(pairs.get("CSN", "").split(":")[0]),
                            "class": unquote(
                                pairs.get("CAVA_CLASS", pairs.get("CLASS", ""))
                            ),
                            "so": unquote(pairs.get("CAVA_SO", pairs.get("SO", ""))),
                            "impact": unquote(
                                pairs.get("CAVA_IMPACT", pairs.get("IMPACT", ""))
                            ),
                            "orig": unquote(pairs.get("CAVA_ORIGHAPLOTYPE", "")),
                            "hap": unquote(pairs.get("CAVA_HAPLOTYPE", "")),
                        }
                    )
            return outrows

    def test_subset_sizes_preserve_exhaustive_cases(self):
        components = ([], ["?"], ["."], ["Ala2Val"], ["Ala2Val", "Gly3Arg"])
        for count, full, split, force, splice in product(
            range(2, 9), components, (False, True), (False, True), (False, True)
        ):
            with self.subTest(count=count, full=full, split=split, force=force, splice=splice):
                sizes = list(haplotype.protein_partition_subset_sizes(count, full, split, force, splice))
                if force or (split and full == ["?"]):
                    self.assertEqual(sizes, [1])
                elif splice or len(full) > 1:
                    self.assertEqual(sizes, list(range(1, count)))
                else:
                    self.assertEqual(sizes, [1])
        for count in (0, 1):
            self.assertEqual(
                list(haplotype.protein_partition_subset_sizes(count, [], True, True, True)),
                [],
            )

    def test_subset_shortcuts_preserve_complete_vcf_output(self):
        tokens = []
        for position in range(7675155, 7675161):
            ref = self.reference.getReference("17", position, position).upper()
            alt = "A" if ref != "A" else "C"
            tokens.append(f"17_{position}_{ref}_{alt}")
        rows = [self._build_row_from_tokens(tokens[:count]) for count in (2, 3, 4, 6)]
        rows += [
            self._build_row_from_tokens([
                "13_32316461_A_C", "13_32316462_T_C", "13_32316463_G_C", "13_32316467_A_C",
            ]),
            self._build_row_from_tokens([
                "13_98463673_A_G", "13_98463692_C_A", "13_98463693_G_A", "13_98463696_G_A",
            ]),
            (
                "chr17", 7675056,
                "chr17_7675054_A_AT;chr17_7675061_TC_T;chr17_7675070_C_CT;chr17_7675076_TG_T",
                "CGCTATCTGAGCAGCGCTCATG", "TCGCTATTGAGCAGCTGCTCAT",
            ),
            (
                "chr17", 7675065,
                "chr17_7675065_A_G;chr17_7675070_C_CT;chr17_7675074_C_T;chr17_7675076_TG_T",
                "AGCAGCGCTCATG", "GGCAGCTGCTTAT",
            ),
        ]
        for split_based, split_adj in product((False, True), repeat=2):
            with self.subTest(split_based=split_based, split_adj=split_adj):
                with redirect_stdout(io.StringIO()):
                    optimized = self._run_rows(rows, split_based, split_adj, return_raw=True)
                    with patch.object(
                        haplotype, "protein_partition_subset_sizes",
                        side_effect=lambda count, *args: range(1, count),
                    ):
                        exhaustive = self._run_rows(rows, split_based, split_adj, return_raw=True)
                self.assertEqual(optimized, exhaustive)

    def test_single_component_shortcut_reduces_annotation_calls(self):
        tokens = []
        for position in range(7675155, 7675161):
            ref = self.reference.getReference("17", position, position).upper()
            alt = "A" if ref != "A" else "C"
            tokens.append(f"17_{position}_{ref}_{alt}")
        row = self._build_row_from_tokens(tokens)
        original_annotate = core.Record.annotate
        with redirect_stdout(io.StringIO()):
            with patch.object(core.Record, "annotate", autospec=True, side_effect=original_annotate) as annotate:
                optimized = self._run_rows([row], return_raw=True)
                optimized_count = annotate.call_count
            with patch.object(
                haplotype, "protein_partition_subset_sizes",
                side_effect=lambda count, *args: range(1, count),
            ), patch.object(core.Record, "annotate", autospec=True, side_effect=original_annotate) as annotate:
                exhaustive = self._run_rows([row], return_raw=True)
                exhaustive_count = annotate.call_count
        self.assertEqual(optimized, exhaustive)
        self.assertEqual(optimized_count, 7)
        self.assertEqual(exhaustive_count, 63)

    def test_overlap_cache_preserves_vcf_and_tsv_output_in_both_lookup_modes(self):
        rows = [
            (
                "chr17", 7675065,
                "chr17_7675065_A_G;chr17_7675070_C_CT;chr17_7675074_C_T;chr17_7675076_TG_T",
                "AGCAGCGCTCATG", "GGCAGCTGCTTAT",
            ),
            self._build_row_from_tokens([
                "13_98463673_A_G", "13_98463692_C_A", "13_98463693_G_A", "13_98463696_G_A",
            ]),
            (
                "17", 43044315,
                "17_43044315_T_A;17_43044317_T_C", "TTT", "ATC",
            ),
        ]
        for load_all, output_format, split_adj in product(
            (False, True), ("VCF", "TSV"), (False, True)
        ):
            with self.subTest(load_all=load_all, output_format=output_format, split_adj=split_adj):
                with redirect_stdout(io.StringIO()):
                    cached = self._run_rows(
                        rows, split_adj=split_adj, return_raw=True,
                        load_all=load_all, output_format=output_format,
                    )
                    with patch.object(
                        Ensembl, "fetch_overlapping_transcripts",
                        Ensembl._fetch_overlapping_transcripts,
                    ):
                        uncached = self._run_rows(
                            rows, split_adj=split_adj, return_raw=True,
                            load_all=load_all, output_format=output_format,
                        )
                self.assertEqual(cached, uncached)

    def test_overlap_cache_reduces_real_annotation_lookup_calls(self):
        row = (
            "chr17", 7675065,
            "chr17_7675065_A_G;chr17_7675070_C_CT;chr17_7675074_C_T;chr17_7675076_TG_T",
            "AGCAGCGCTCATG", "GGCAGCTGCTTAT",
        )
        original = Ensembl._fetch_overlapping_transcripts
        for load_all in (False, True):
            with self.subTest(load_all=load_all), redirect_stdout(io.StringIO()):
                with patch.object(Ensembl, "_fetch_overlapping_transcripts", autospec=True, side_effect=original) as fetch:
                    cached = self._run_rows([row], split_adj=True, return_raw=True, load_all=load_all)
                    cached_count = fetch.call_count
                with patch.object(Ensembl, "fetch_overlapping_transcripts", autospec=True, side_effect=original) as fetch:
                    uncached = self._run_rows([row], split_adj=True, return_raw=True, load_all=load_all)
                    uncached_count = fetch.call_count
                self.assertEqual(cached, uncached)
                self.assertEqual(cached_count, 6)
                self.assertEqual(uncached_count, 64)

    def test_four_variant_splice_impact_is_not_worse_than_components(self):
        donor_side_tokens = [
            "13_32316525_C_A",
            "13_32316526_A_C",
            "13_32316529_T_A",
            "13_32316533_G_A",
        ]
        acceptor_side_tokens = [
            "13_98463673_A_G",
            "13_98463692_C_A",
            "13_98463693_G_A",
            "13_98463696_G_A",
        ]

        rows = [
            self._build_row_from_tokens(donor_side_tokens),
            self._build_row_from_tokens(acceptor_side_tokens),
        ]
        outrows = self._run_rows(rows, split_based=False, split_adj=False)
        by_orig = defaultdict(list)
        for r in outrows:
            if r["orig"]:
                by_orig[r["orig"]].append(r)

        for tokens in (donor_side_tokens, acceptor_side_tokens):
            orig = ";".join(tokens)
            emitted = by_orig.get(orig, [])
            self.assertGreaterEqual(len(emitted), 2, msg=orig)

            singleton_rows = [
                (tok.split("_", 3)[0], int(tok.split("_", 3)[1]), tok, tok.split("_", 3)[2], tok.split("_", 3)[3])
                for tok in tokens
            ]
            singleton_out = self._run_rows(
                singleton_rows, split_based=False, split_adj=False
            )
            singleton_severity = [
                self._splice_severity(r["class"], r["so"]) for r in singleton_out
            ]
            max_singleton_severity = max(singleton_severity)

            singleton_splice_impacts = []
            for r, sev in zip(singleton_out, singleton_severity):
                if sev > 0 and r["impact"].isdigit():
                    singleton_splice_impacts.append(int(r["impact"]))
            min_singleton_splice_impact = (
                min(singleton_splice_impacts) if singleton_splice_impacts else None
            )

            for rec in emitted:
                rec_sev = self._splice_severity(rec["class"], rec["so"])
                self.assertLessEqual(rec_sev, max_singleton_severity, msg=orig)
                if (
                    rec_sev > 0
                    and min_singleton_splice_impact is not None
                    and rec["impact"].isdigit()
                ):
                    self.assertGreaterEqual(
                        int(rec["impact"]), min_singleton_splice_impact, msg=orig
                    )

        acceptor_orig = ";".join(acceptor_side_tokens)
        self.assertFalse(
            any(r["hap"] == acceptor_orig for r in by_orig[acceptor_orig]),
            msg=acceptor_orig,
        )

    def test_four_variant_first_coding_base_combination_splits(self):
        tokens = [
            "13_32316461_A_C",
            "13_32316462_T_C",
            "13_32316463_G_C",
            "13_32316467_A_C",
        ]
        row = self._build_row_from_tokens(tokens)
        outrows = self._run_rows([row], split_based=True, split_adj=False)

        emitted = [r for r in outrows if r["orig"] == row[2]]
        self.assertGreaterEqual(len(emitted), 2)

        split_csns = [r["csn"] for r in emitted if r["hap"] != row[2]]
        self.assertIn("c.1A>C_p.Met1?", split_csns)
        self.assertIn("c.2T>C_p.Met1?", split_csns)
        self.assertIn("c.3G>C_p.Met1?", split_csns)

    def test_four_variant_mixed_proximity_can_partition_into_near_and_far_blocks(self):
        row = (
            "chr17",
            7675065,
            "chr17_7675065_A_G;chr17_7675070_C_CT;chr17_7675074_C_T;chr17_7675076_TG_T",
            "AGCAGCGCTCATG",
            "GGCAGCTGCTTAT",
        )
        outrows = self._run_rows([row], split_based=True, split_adj=False)

        emitted = [r for r in outrows if r["orig"] == row[2] and r["hap"] != row[2]]
        self.assertGreaterEqual(len(emitted), 2)

        subset_sizes = [len([x for x in r["hap"].split(";") if x]) for r in emitted]
        self.assertIn(1, subset_sizes)
        self.assertIn(3, subset_sizes)

    def test_parse_mode_splits_splice_spanning_haplotypes_without_atomic_splice_support(
        self,
    ):
        rows = [
            (
                "13",
                98463673,
                "chr13_98463673_A_G;chr13_98463693_G_A",
                "ACTCAGGGGCCGACTTACGCG",
                "GCTCAGGGGCCGACTTACGCA",
            ),
            (
                "13",
                98463673,
                "chr13_98463673_A_G;chr13_98463696_G_A",
                "ACTCAGGGGCCGACTTACGCGTCG",
                "GCTCAGGGGCCGACTTACGCGTCA",
            ),
            (
                "13",
                98463673,
                "chr13_98463673_A_G;chr13_98463693_G_A;chr13_98463696_G_A",
                "ACTCAGGGGCCGACTTACGCGTCG",
                "GCTCAGGGGCCGACTTACGCATCA",
            ),
        ]

        outrows = self._run_rows(rows, split_based=False, split_adj=False)
        by_orig = defaultdict(list)
        for r in outrows:
            if r["orig"]:
                by_orig[r["orig"]].append(r)

        for _, _, row_id, _, _ in rows:
            emitted = by_orig.get(row_id, [])
            self.assertGreaterEqual(len(emitted), 2, msg=row_id)
            self.assertFalse(any(r["hap"] == row_id for r in emitted), msg=row_id)

            orig_n = len([x for x in row_id.split(";") if x])
            for rec in emitted:
                split_n = len([x for x in rec["hap"].split(";") if x])
                self.assertLess(split_n, orig_n, msg=row_id)
                self.assertNotIn("ESS", rec["class"], msg=row_id)
                self.assertNotIn("splice_acceptor_variant", rec["so"], msg=row_id)

    def test_adjacent_and_slightly_separated_start_codon_haplotypes_split_to_singletons(
        self,
    ):
        # These atomic variants are already covered as single-variant tests in test_end2end.py.
        atoms = [
            ("13", 32316461, "13_32316461_A_C", "A", "C", "c.1A>C_p.Met1?"),
            ("13", 32316462, "13_32316462_T_C", "T", "C", "c.2T>C_p.Met1?"),
            ("13", 32316463, "13_32316463_G_C", "G", "C", "c.3G>C_p.Met1?"),
        ]

        rows = [a[:5] for a in atoms] + [c["row"] for c in START_CODON_MULTI_CASES]
        outrows = self._run_rows(rows)

        by_orig = defaultdict(list)
        for r in outrows:
            if r["orig"]:
                by_orig[r["orig"]].append(r)

        for case in START_CODON_MULTI_CASES:
            records = by_orig[case["row"][2]]
            split_csns = [r["csn"] for r in records if r["hap"] != case["row"][2]]
            self.assertEqual(
                sorted(case["expected_split"]), sorted(split_csns), msg=case["name"]
            )

    def test_all_two_serine_snv_pairs_split_to_independent_expected_arrays(self):
        # Known consecutive Ser codons from existing single-variant region:
        # chr17:43092111-43092116 is ACTACT on genome (+), coding Ser1140/Ser1139 on transcript (-).
        chrom = "17"
        start = 43092111
        ref_seq = "ACTACT"

        # Build all single SNVs across the 6 bp span.
        single_rows = []
        for i, ref_base in enumerate(ref_seq):
            pos = start + i
            for alt_base in "ACGT":
                if alt_base == ref_base:
                    continue
                token = f"{chrom}_{pos}_{ref_base}_{alt_base}"
                single_rows.append((chrom, pos, token, ref_base, alt_base))

        # First pass: annotate singles to confirm the explicit expected array below matches
        # the already-existing single-variant tests for each atomic SNV.
        single_out = self._run_rows(single_rows)
        for r in single_out:
            if not r["orig"]:
                self.assertEqual(
                    SERINE_SINGLE_EXPECTATIONS[r["id"]], r["csn"], msg=r["id"]
                )

        # Group atoms by codon block (one SNV from each consecutive Ser codon).
        left_block = [r for r in single_rows if 43092111 <= r[1] <= 43092113]
        right_block = [r for r in single_rows if 43092114 <= r[1] <= 43092116]

        def build_haplotype_row(a, b):
            chrom_a, pos_a, token_a, _, alt_a = a
            chrom_b, pos_b, token_b, _, alt_b = b
            self.assertEqual(chrom_a, chrom_b)
            lo = min(pos_a, pos_b)
            hi = max(pos_a, pos_b)
            ref_span = list(ref_seq[lo - start : hi - start + 1])
            alt_span = ref_span[:]
            alt_span[pos_a - lo] = alt_a
            alt_span[pos_b - lo] = alt_b
            if pos_a < pos_b:
                row_id = token_a + ";" + token_b
                expected = [
                    SERINE_SINGLE_EXPECTATIONS[token_a],
                    SERINE_SINGLE_EXPECTATIONS[token_b],
                ]
            else:
                row_id = token_b + ";" + token_a
                expected = [
                    SERINE_SINGLE_EXPECTATIONS[token_b],
                    SERINE_SINGLE_EXPECTATIONS[token_a],
                ]
            return (chrom_a, lo, row_id, "".join(ref_span), "".join(alt_span)), expected

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
            by_orig[r["orig"]].append(r)

        synonymous_seen = 0
        for row_id, expected in expected_arrays:
            records = by_orig[row_id]
            split_records = [r for r in records if r["hap"] != row_id]
            split_csns = [r["csn"] for r in split_records]
            self.assertEqual(sorted(expected), sorted(split_csns), msg=row_id)
            synonymous_seen += sum(
                1 for e in expected if "_p.Ser" in e and e.endswith("=")
            )

        self.assertGreater(synonymous_seen, 0)

    def test_mixed_splice_region_and_coding_haplotype_suppresses_canonical(self):
        row = (
            "chr17",
            7675056,
            "chr17_7675054_A_AT;chr17_7675061_TC_T;chr17_7675070_C_CT;chr17_7675076_TG_T",
            "CGCTATCTGAGCAGCGCTCATG",
            "TCGCTATTGAGCAGCTGCTCAT",
        )
        outrows = self._run_rows([row], split_based=True, split_adj=True)

        emitted = [r for r in outrows if r["orig"] == row[2]]
        self.assertGreaterEqual(len(emitted), 2)
        self.assertFalse(any(r["hap"] == row[2] for r in emitted))
        self.assertTrue(any("splice_region_variant" in r["so"] for r in emitted))
        self.assertTrue(any("frameshift_variant" in r["so"] for r in emitted))

    def test_utr_variants_one_bp_apart_force_split(self):
        row = (
            "17",
            43044315,
            "17_43044315_T_A;17_43044317_T_C",
            "TTT",
            "ATC",
        )
        outrows = self._run_rows([row], split_based=False, split_adj=False)

        emitted = [r for r in outrows if r["orig"] == row[2]]
        self.assertGreaterEqual(len(emitted), 2)
        self.assertFalse(any(r["hap"] == row[2] for r in emitted))
        self.assertTrue(all(r["class"] == "3PU" for r in emitted))

    def test_adjacent_utr_variants_do_not_force_split(self):
        row = (
            "17",
            43044315,
            "17_43044315_T_A;17_43044316_T_C",
            "TT",
            "AC",
        )
        outrows = self._run_rows([row], split_based=False, split_adj=False)

        emitted = [r for r in outrows if r["orig"] == row[2]]
        self.assertEqual(1, len(emitted))
        self.assertEqual(row[2], emitted[0]["hap"])
        self.assertEqual("3PU", emitted[0]["class"])


if __name__ == "__main__":
    unittest.main()
