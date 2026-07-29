import json
import os
import tempfile
import unittest
from collections import defaultdict
from urllib.parse import unquote

from cava.utils import core
from cava.utils import haplotype
from cava.utils import main
from cava.utils.data import Reference


class TestHaplotypePromptJsonConformance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        cls.workspace_root = os.path.dirname(cls.repo_root)
        cls.fixture_json = os.path.join(
            cls.workspace_root,
            "Prompt_for_Cis",
            "hgvs_cis_haplotype_GRCh38_final",
            "hgvs_cis_unit_tests_grch38_mane.json",
        )

        with open(cls.fixture_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        cls.all_rows = []
        for section in (
            "two_variant_tests",
            "three_variant_tests",
            "four_variant_tests",
        ):
            cls.all_rows.extend(payload.get(section, []))

        cls.rows = [
            r
            for r in cls.all_rows
            if r.get("VCF_SEQUENCE_STATUS") == "VERIFIED_GRCH38_PLUS_STRAND"
            and r.get("test_id") != "2V-031"
        ]

    def _run_single_job(self, rows, split_based=True, split_adj=True):
        with tempfile.TemporaryDirectory() as td:
            cfg = os.path.join(td, "cfg.txt")
            inp = os.path.join(td, "in.vcf")
            outprefix = os.path.join(td, "out")

            with open(cfg, "w", encoding="utf-8") as f:
                f.write("@inputformat=VCF\n")
                f.write("@outputformat=VCF\n")
                f.write(
                    "@reference="
                    + os.path.join(self.repo_root, "cava", "data", "tmp.GRCh38.fa")
                    + "\n"
                )
                f.write(
                    "@ensembl="
                    + os.path.join(
                        self.repo_root,
                        "cava",
                        "data",
                        "MANE.GRCh38.v1.1.refseq_genomic.db.gz",
                    )
                    + "\n"
                )
                f.write("@dbsnp=.\n")
                f.write("@logfile=FALSE\n")
                f.write("@prefix=FALSE\n")
                f.write("@chrom=.\n")

            reference = Reference(core.Options(cfg))

            def materialize_row(row):
                pos = row["VCFPOS"]
                ref = row["VCFREF"]
                alt = row["VCFALT"]
                if not str(pos).strip() or not str(ref).strip() or not str(alt).strip():
                    atoms = [
                        haplotype.parse_atomic_token(tok)
                        for tok in row["VCFID"].split(";")
                        if tok.strip()
                    ]
                    pos, ref, alt, _ = haplotype.build_subset_vcf_fields(
                        reference, row["VCFCHROM"], atoms
                    )
                return row["VCFCHROM"], int(pos), row["VCFID"], ref, alt

            unique = {}
            for r in rows:
                key = materialize_row(r)
                unique[key] = r

            with open(inp, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for chrom, pos, vid, ref, alt in unique.keys():
                    f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")

            options = core.Options(cfg)
            options.args["parseHaplotype"] = True
            options.args["splitBasedOnProtein"] = bool(split_based)
            options.args["splitadjacentprotein"] = bool(split_adj)

            with open(outprefix + ".vcf", "w", encoding="utf-8") as out:
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
                1, options, copts, 3, "", set(), set(), set(), impactdir, len(unique)
            )
            job.run()

            outrows = []
            with open(outprefix + ".vcf", "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    outrows.append(line.strip().split("\t"))
            return outrows

    def test_all_prompt_rows_have_expected_csn_output(self):
        outrows = self._run_single_job(self.rows, split_based=True, split_adj=True)
        by_orig = defaultdict(set)

        for cols in outrows:
            info = cols[7]
            pairs = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    pairs[k] = v
            csn = unquote(pairs.get("CSN", "").split(":")[0])
            protein = csn.split("_p.", 1)[1] if "_p." in csn else csn
            orig = unquote(pairs.get("CAVA_ORIGHAPLOTYPE", ""))
            if orig:
                for comp in haplotype.protein_components_from_expected_p_hgvs(
                    "p." + protein
                ):
                    by_orig[orig].add(comp)

        missing = []
        for r in self.rows:
            observed = by_orig.get(r["VCFID"], set())
            expected_components = haplotype.protein_components_from_expected_p_hgvs(
                r["expected_p_hgvs"]
            )
            for comp in expected_components:
                if comp not in observed:
                    missing.append((r["test_id"], comp))

        self.assertEqual([], missing)

    def test_splice_p_unknown_still_splits_nearby(self):
        row = next(
            r
            for r in self.all_rows
            if r["test_id"] == "2V-013"
        )
        outrows = self._run_single_job([row], split_based=True, split_adj=False)

        self.assertGreaterEqual(len(outrows), 2)
        orig = None
        saw_subset = False
        for cols in outrows:
            info = cols[7]
            pairs = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    pairs[k] = v
            if orig is None:
                orig = pairs.get("CAVA_ORIGHAPLOTYPE")
            if pairs.get("CAVA_HAPLOTYPE") != orig:
                saw_subset = True
        self.assertTrue(saw_subset)


if __name__ == "__main__":
    unittest.main()
