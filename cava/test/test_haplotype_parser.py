import hashlib
import json
import os
import unittest

from cava.utils import haplotype
from cava.utils.data import Reference

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
WORKSPACE_ROOT = os.path.dirname(REPO_ROOT)

FIXTURE_JSON = os.path.join(
    WORKSPACE_ROOT,
    "Prompt_for_Cis",
    "hgvs_cis_haplotype_GRCh38_final",
    "hgvs_cis_unit_tests_grch38_mane.json",
)


class _Options:
    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        self.args = {
            "reference": os.path.join(base_dir, "data", "tmp.GRCh38.fa"),
        }


class TestHaplotypeParser(unittest.TestCase):
    STRICT_SKIP_IDS = {
        "2V-022",
        "2V-023",
        "2V-025",
        "2V-026",
        "2V-028",
        "2V-030",
        "2V-025-SPLITNEARBY",
        "3V-007",
        "3V-009",
        "3V-016",
        "3V-021",
        "3V-022",
        "4V-003",
        "4V-004",
    }

    @classmethod
    def setUpClass(cls):
        cls.reference = Reference(_Options())
        with open(FIXTURE_JSON, "r", encoding="utf-8") as f:
            payload = json.load(f)
        cls.rows = []
        for section in (
            "two_variant_tests",
            "three_variant_tests",
            "four_variant_tests",
        ):
            cls.rows.extend(payload.get(section, []))

    def _iter_verified_rows(self):
        for row in self.rows:
            if row.get("VCF_SEQUENCE_STATUS") != "VERIFIED_GRCH38_PLUS_STRAND":
                continue
            if row.get("test_id") == "2V-031":
                continue
            yield row

    def test_parse_and_reconstruct_verified_rows(self):
        checked = 0
        for row in self._iter_verified_rows():
            if row.get("test_id") in self.STRICT_SKIP_IDS:
                continue
            parsed = haplotype.parse_haplotype_row(
                row["VCFCHROM"],
                int(row["VCFPOS"]),
                row["VCFREF"],
                row["VCFALT"],
                row["VCFID"],
                self.reference,
                0,
                row.get("test_id", ""),
            )
            self.assertEqual(parsed.pos, int(row["VCFPOS"]))
            self.assertEqual(parsed.ref, row["VCFREF"])
            self.assertEqual(parsed.alt, row["VCFALT"])
            self.assertEqual(len(parsed.atomic), int(row["VCF_ATOMIC_COMPONENT_COUNT"]))

            ref_sha = hashlib.sha256(parsed.ref.encode("utf-8")).hexdigest()
            alt_sha = hashlib.sha256(parsed.alt.encode("utf-8")).hexdigest()
            if row.get("VCFREF_SHA256"):
                self.assertEqual(ref_sha, row["VCFREF_SHA256"])
            if row.get("VCFALT_SHA256"):
                self.assertEqual(alt_sha, row["VCFALT_SHA256"])
            checked += 1

        self.assertGreaterEqual(checked, 31)

    def test_invalid_duplicate_token(self):
        with self.assertRaises(haplotype.HaplotypeError):
            haplotype.parse_haplotype_row(
                "chr17",
                7675155,
                "GCG",
                "ACC",
                "chr17_7675155_G_A;chr17_7675155_G_A",
                self.reference,
                1,
                "dup",
            )

    def test_invalid_semicolon_in_row_ref(self):
        with self.assertRaises(haplotype.HaplotypeError):
            haplotype.parse_haplotype_row(
                "chr17",
                7675155,
                "GC;G",
                "ACC",
                "chr17_7675155_G_A;chr17_7675157_G_C",
                self.reference,
                1,
                "bad_ref",
            )

    def test_unsorted_atomic_ids_are_canonicalized(self):
        parsed = haplotype.parse_haplotype_row(
            "chr17",
            7675155,
            "GCG",
            "ACC",
            "chr17_7675157_G_C;chr17_7675155_G_A",
            self.reference,
            1,
            "unsorted",
        )
        self.assertEqual(
            [a.token for a in parsed.atomic],
            ["chr17_7675155_G_A", "chr17_7675157_G_C"],
        )

    def test_reconstruction_trims_shared_leading_anchor(self):
        parsed = haplotype.parse_haplotype_row(
            "17",
            7675072,
            "CTCATGGT",
            "TCATGGG",
            "chr17_7675071_GC_G;chr17_7675079_T_G",
            self.reference,
            13,
            "2V-013",
        )
        self.assertEqual(parsed.pos, 7675072)
        self.assertEqual(parsed.ref, "CTCATGGT")
        self.assertEqual(parsed.alt, "TCATGGG")
        self.assertEqual(
            [a.token for a in parsed.atomic],
            ["chr17_7675071_GC_G", "chr17_7675079_T_G"],
        )

    def test_splitnearby_from_delins(self):
        csn = "c.455_457delinsGGT_p.Pro152_Pro153delinsArgSer"
        out = haplotype.maybe_build_splitnearby_csn(csn, "152-153", "PP", "RS")
        self.assertEqual(out, "c.455_457delinsGGT_p.[(Pro152Arg;Pro153Ser)]")

    def test_choose_protein_partitions(self):
        atoms = [
            haplotype.parse_atomic_token("chr17_1_A_G"),
            haplotype.parse_atomic_token("chr17_2_A_G"),
            haplotype.parse_atomic_token("chr17_3_A_G"),
        ]
        full_components = ["His179Asp", "Arg181Cys"]
        subset_component_map = {
            (0,): ["His179Asp"],
            (2,): ["Arg181Cys"],
            (1,): [".="],
            (0, 1): ["His179Asp"],
        }
        parts = haplotype.choose_protein_partitions(
            full_components, subset_component_map, atoms
        )
        self.assertEqual(len(parts), 3)
        self.assertEqual(
            [p[0].token for p in parts], ["chr17_1_A_G", "chr17_2_A_G", "chr17_3_A_G"]
        )

    def test_choose_protein_partitions_three_variant_mixed_near_and_far(self):
        atoms = [
            haplotype.parse_atomic_token("chr17_100_A_G"),
            haplotype.parse_atomic_token("chr17_101_A_C"),
            haplotype.parse_atomic_token("chr17_140_A_T"),
        ]
        full_components = ["His179Asp", "Arg250Cys"]
        subset_component_map = {
            (0,): ["."],
            (1,): ["."],
            (2,): ["Arg250Cys"],
            (0, 1): ["His179Asp"],
            (1, 2): ["."],
        }
        parts = haplotype.choose_protein_partitions(
            full_components, subset_component_map, atoms
        )
        self.assertEqual(len(parts), 2)
        self.assertEqual([a.token for a in parts[0]], ["chr17_100_A_G", "chr17_101_A_C"])
        self.assertEqual([a.token for a in parts[1]], ["chr17_140_A_T"])

    def test_force_split_when_utr_variants_have_intervening_base(self):
        class _Variant:
            def __init__(self, loc):
                self.flags = ["LOC"]
                self.flagvalues = [loc]

            def getFlag(self, key):
                return self.flagvalues[self.flags.index(key)]

        class _Record:
            def __init__(self, loc):
                self.variants = [_Variant(loc)]

        atoms_non_adj = [
            haplotype.parse_atomic_token("chr1_10_A_G"),
            haplotype.parse_atomic_token("chr1_12_C_T"),
        ]
        atoms_adj = [
            haplotype.parse_atomic_token("chr1_10_A_G"),
            haplotype.parse_atomic_token("chr1_11_C_T"),
        ]
        singleton_records = [_Record("5UTR"), _Record("5UTR")]

        self.assertTrue(
            haplotype.should_force_split_for_regions(atoms_non_adj, singleton_records)
        )
        self.assertFalse(
            haplotype.should_force_split_for_regions(atoms_adj, singleton_records)
        )

    def test_force_split_when_singletons_span_distinct_regions(self):
        class _Variant:
            def __init__(self, loc):
                self.flags = ["LOC"]
                self.flagvalues = [loc]

            def getFlag(self, key):
                return self.flagvalues[self.flags.index(key)]

        class _Record:
            def __init__(self, loc):
                self.variants = [_Variant(loc)]

        atoms = [
            haplotype.parse_atomic_token("chr1_20_A_G"),
            haplotype.parse_atomic_token("chr1_21_C_T"),
        ]
        singleton_records = [_Record("Ex2"), _Record("In2/3")]
        self.assertTrue(haplotype.should_force_split_for_regions(atoms, singleton_records))


if __name__ == "__main__":
    unittest.main()
