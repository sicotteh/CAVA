import hashlib
import json
import os
import unittest
from itertools import product
from random import Random

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

    def test_minimal_vcf_matches_original_algorithm_exhaustively(self):
        def original_minimal(pos, ref, alt):
            while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
                ref = ref[:-1]
                alt = alt[:-1]
            while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                pos += 1
            return pos, ref, alt

        alleles = [
            "".join(bases)
            for length in range(7)
            for bases in product("AC", repeat=length)
        ]
        for ref, alt in product(alleles, repeat=2):
            self.assertEqual(
                haplotype._minimal_vcf(7, ref, alt),
                original_minimal(7, ref, alt),
                (ref, alt),
            )
        for ref, alt in (
            ("A" * 50000 + "C", "A" * 50000 + "G"),
            ("C" + "A" * 50000, "G" + "A" * 50000),
            ("A" * 1000 + "C" + "T" * 1000, "A" * 1000 + "G" + "T" * 1000),
        ):
            self.assertEqual(
                haplotype._minimal_vcf(7, ref, alt), original_minimal(7, ref, alt)
            )

    def test_reconstruction_does_not_refetch_validated_trimmed_reference(self):
        class CountingReference:
            def __init__(self):
                self.calls = []

            def getReference(self, chrom, start, end):
                self.calls.append((chrom, start, end))
                return "AACCGGTT"[start - 1:end]

        reference = CountingReference()
        atoms = [
            haplotype.parse_atomic_token("1_1_AA_AC"),
            haplotype.parse_atomic_token("1_4_CG_C"),
            haplotype.parse_atomic_token("1_7_T_TA"),
        ]
        self.assertEqual(
            haplotype._reconstruct_from_atoms(reference, "1", atoms),
            (2, "ACCGGT", "CCCGTA"),
        )
        self.assertEqual(len(reference.calls), len(atoms) + 1)
        self.assertEqual(reference.calls[-1], ("1", 1, 7))

    def test_reconstruction_still_validates_full_atomic_reference(self):
        class Reference:
            def getReference(self, chrom, start, end):
                return "AACCGGTT"[start - 1:end]

        for token in ("1_1_CA_CC", "1_3_CA_GA", "1_4_TG_T"):
            with self.subTest(token=token):
                with self.assertRaisesRegex(haplotype.HaplotypeError, "Atomic REF mismatch"):
                    haplotype._reconstruct_from_atoms(
                        Reference(), "1", [haplotype.parse_atomic_token(token)]
                    )

    def test_reconstruction_matches_original_overlap_and_edit_application(self):
        genome = "ACAC"

        class Reference:
            def getReference(self, chrom, start, end):
                return genome[start - 1:end]

        def original_reconstruction(atoms):
            edits = []
            for atom in atoms:
                start, end, _, replacement = haplotype._trim_atomic_edit(atom)
                edits.append((start - 1, end, replacement, atom.token))
            spans = [edit for edit in edits if edit[1] > edit[0]]
            for index, left in enumerate(spans):
                for right in spans[index + 1:]:
                    if max(left[0], right[0]) < min(left[1], right[1]):
                        raise haplotype.HaplotypeError(
                            f"Overlapping or unsorted atomic edits: {left[3]} and {right[3]}"
                        )
            start = min(min(edit[0] for edit in edits), min(atom.pos - 1 for atom in atoms))
            end = max(max(edit[1] for edit in edits), max(atom.pos - 1 + len(atom.ref) for atom in atoms))
            ref = genome[start:end]
            alt = ref
            for edit_start, edit_end, replacement, _ in sorted(
                edits, key=lambda edit: (edit[0], edit[1]), reverse=True
            ):
                alt = alt[:edit_start - start] + replacement + alt[edit_end - start:]
            return haplotype._minimal_vcf(start + 1, ref, alt)

        alternatives = [
            "".join(bases)
            for length in (1, 2)
            for bases in product("AC", repeat=length)
        ]
        candidates = []
        for position in range(1, len(genome) + 1):
            for ref_length in (1, 2):
                if position - 1 + ref_length > len(genome):
                    continue
                ref = genome[position - 1:position - 1 + ref_length]
                for alt in alternatives:
                    candidates.append(haplotype.parse_atomic_token(f"1_{position}_{ref}_{alt}"))
        reference = Reference()
        checked = 0
        for size in (1, 2, 3):
            for atoms in product(candidates, repeat=size):
                try:
                    expected = original_reconstruction(atoms)
                except haplotype.HaplotypeError as error:
                    with self.assertRaises(haplotype.HaplotypeError) as actual:
                        haplotype._reconstruct_from_atoms(reference, "1", atoms)
                    self.assertEqual(str(actual.exception), str(error))
                else:
                    self.assertEqual(
                        haplotype._reconstruct_from_atoms(reference, "1", atoms),
                        expected,
                        tuple(atom.token for atom in atoms),
                    )
                checked += 1
        self.assertEqual(checked, 42 + 42 ** 2 + 42 ** 3)
        random = Random(417)
        for _ in range(2000):
            genome = "".join(random.choices("ACGT", k=80))
            atoms = []
            for _ in range(random.randint(1, 12)):
                position = random.randint(1, len(genome) - 6)
                ref_length = random.randint(1, 6)
                ref = genome[position - 1:position - 1 + ref_length]
                alt = "".join(random.choices("ACGT", k=random.randint(1, 8)))
                atoms.append(haplotype.parse_atomic_token(f"1_{position}_{ref}_{alt}"))
            try:
                expected = original_reconstruction(atoms)
            except haplotype.HaplotypeError as error:
                with self.assertRaises(haplotype.HaplotypeError) as actual:
                    haplotype._reconstruct_from_atoms(reference, "1", atoms)
                self.assertEqual(str(actual.exception), str(error))
            else:
                self.assertEqual(
                    haplotype._reconstruct_from_atoms(reference, "1", atoms),
                    expected,
                    tuple(atom.token for atom in atoms),
                )

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
