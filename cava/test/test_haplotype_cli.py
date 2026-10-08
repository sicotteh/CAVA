import io
import os
import re
import tempfile
import unittest
from contextlib import redirect_stdout
from itertools import product

from cava.utils import core
from cava.utils import haplotype
from cava.utils import main
from cava.utils.data import Reference


class _Options:
    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        self.args = {
            "reference": os.path.join(base_dir, "data", "tmp.GRCh38.fa"),
        }


class TestHaplotypeCLI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = Reference(_Options())
        cls.repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    def test_header_has_haplotype_tags(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        out = tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8")
        try:
            core.writeHeader(options, "", out, False, "2.0.15")
            out.close()
            with open(out.name, "r", encoding="utf-8") as f:
                text = f.read()
            self.assertIn("##INFO=<ID=CAVA_ORIGHAPLOTYPE", text)
            self.assertIn("##INFO=<ID=CAVA_HAPLOTYPE", text)
        finally:
            try:
                os.unlink(out.name)
            except OSError:
                pass

    def test_header_defines_csnalt_as_single_string(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        for prefix in (False, True):
            with self.subTest(prefix=prefix):
                options.args["prefix"] = prefix
                out = io.StringIO()
                core.writeHeader(options, "", out, False, "2.0.15")
                key = "CAVA_CSNALT" if prefix else "CSNALT"
                self.assertIn(
                    "##INFO=<ID=" + key + ",Number=1,Type=String,",
                    out.getvalue(),
                )

    def test_record_output_encodes_annotation_semicolons(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        annotation = "c.477_486delinsACGT_p.[(Met159Met;Ala160Ala;Val161Val;Thr162Ala)]"
        for prefix in (False, True):
            with self.subTest(prefix=prefix):
                options.args["prefix"] = prefix
                record = core.Record(
                    "chr17\t7675155\t.\tG\tA\t.\tPASS\tCSNALT=old;CAVA_CSNALT=old;DP=7\n",
                    options, None, self.reference,
                )
                for key in ("CSN", "CSNALT", "ALTANN"):
                    record.variants[0].addFlag(key, annotation)
                record.variants[0].addFlag("TRANSCRIPT", "NM_TEST")
                record.variants[0].addFlag("GENE", "TEST")
                options.transcript2protein = {"NM_TEST": "NP_TEST"}
                record.haplotype_hgvsc_override = "c.[477G>A;486T>C]"
                record.variants[0].addFlag("CUSTOM", "already%3Bencoded;raw")
                out = io.StringIO()
                record.output("VCF", out, options, set(), set(), set(), False)
                fields = out.getvalue().strip().split("\t")[7].split(";")
                self.assertTrue(all("=" in field for field in fields))
                pairs = dict(field.split("=", 1) for field in fields)
                tag_prefix = "CAVA_" if prefix else ""
                for key in ("CSN", "CSNALT", "ALTANN"):
                    self.assertEqual(
                        pairs[tag_prefix + key], annotation.replace(";", "%3B")
                    )
                self.assertEqual(pairs[tag_prefix + "CUSTOM"], "already%3Bencoded%3Braw")
                self.assertEqual(
                    pairs[tag_prefix + "HGVSc"], "c.[477G>A%3B486T>C]"
                )
                self.assertEqual(
                    pairs[tag_prefix + "HGVSp"],
                    "NP_TEST:p.([(Met159Met%3BAla160Ala%3BVal161Val%3BThr162Ala)])",
                )
                self.assertNotIn("old", out.getvalue())
                self.assertEqual(pairs["DP"], "7")
                self.assertEqual(record.variants[0].getFlag("CSNALT"), annotation)

    def test_vcf_headers_and_tsv_fields_match_emitted_annotations(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        options.args["ensembl"] = "test.db.gz"
        options.args["dbsnp"] = "test.dbsnp.gz"
        options.transcript2protein = {"NM_ONE": "NP_ONE", "NM_TWO": "NP_TWO"}
        annotations = {
            "TRANSCRIPT": "NM_ONE:NM_TWO", "GENE": "ONE:TWO",
            "GENEID": "1:2", "TRINFO": "info1:info2", "LOC": "Ex1:Ex2",
            "CSN": "c.1A>G_p.Met1Val:c.2C>T_p.Ala2Val",
            "PROTPOS": "1:2", "PROTREF": "M:A", "PROTALT": "V:V",
            "CSNALT": "c.1A>G_p.[(Met1Val;Ala2Ala)]:c.2C>T_p.Ala2Val",
            "CLASS": "NSY:NSY", "SO": "missense_variant:missense_variant",
            "IMPACT": "2:2", "ALTANN": "alt1:alt2", "ALTCLASS": "NSY:NSY",
            "ALTSO": "missense_variant:missense_variant", "ALTFLAG": "None:None",
            "DBSNP": "rs123", "HGVSg": "NC_000017.11:g.7675155G>A",
            "CAVA_ORIGHAPLOTYPE": "atom1%3Batom2", "CAVA_HAPLOTYPE": "atom1",
        }
        combinations = product((False, True), ("CLASS", "SO", "BOTH"), (False, True), (False, True))
        for prefix, ontology, givealt, givealtflag in combinations:
            with self.subTest(prefix=prefix, ontology=ontology, givealt=givealt, givealtflag=givealtflag):
                options.args.update(prefix=prefix, ontology=ontology, givealt=givealt, givealtflag=givealtflag)
                record = core.Record(
                    "chr17\t7675155\t.\tG\tA,C\t.\tPASS\t.\n",
                    options, None, self.reference,
                )
                for variant in record.variants:
                    for key, value in reversed(list(annotations.items())):
                        if key in {"CLASS", "ALTCLASS"} and ontology == "SO":
                            continue
                        if key in {"SO", "ALTSO"} and ontology == "CLASS":
                            continue
                        if key in {"ALTANN", "ALTCLASS", "ALTSO"} and not givealt:
                            continue
                        if key == "ALTFLAG" and givealt and not givealtflag:
                            continue
                        variant.addFlag(key, value)

                options.args["outputformat"] = "VCF"
                header = io.StringIO()
                core.writeHeader(options, "", header, False, "2.0.15")
                defined = set(re.findall(r"##INFO=<ID=([^,]+),", header.getvalue()))
                outvcf = io.StringIO()
                record.output("VCF", outvcf, options, set(), set(), set(), False)
                info = dict(field.split("=", 1) for field in outvcf.getvalue().strip().split("\t")[7].split(";"))
                self.assertFalse(set(info) - defined)

                options.args["outputformat"] = "TSV"
                header = io.StringIO()
                core.writeHeader(options, "", header, False, "2.0.15")
                columns = header.getvalue().strip().split("\t")
                outtsv = io.StringIO()
                record.output("TSV", outtsv, options, set(), set(), set(), False)
                rows = [line.split("\t") for line in outtsv.getvalue().splitlines()]
                self.assertEqual(len(rows), 4)
                for row in rows:
                    self.assertEqual(len(row), len(columns))
                column_keys = {"HGVSG": "HGVSg", "HGVSC": "HGVSc", "HGVSP": "HGVSp"}
                for key, value in info.items():
                    raw_key = key[5:] if prefix and key.startswith("CAVA_") and key not in {"CAVA_ORIGHAPLOTYPE", "CAVA_HAPLOTYPE"} else key
                    column = next((name for name in columns if column_keys.get(name, name) == raw_key), None)
                    self.assertIsNotNone(column, key)
                    index = columns.index(column)
                    if raw_key in {"TYPE", "DBSNP", "HGVSg", "CAVA_ORIGHAPLOTYPE", "CAVA_HAPLOTYPE"}:
                        expected = value.split(",")[0]
                        self.assertTrue(all(row[index].replace(";", "%3B") == expected for row in rows), key)
                    else:
                        actual = ",".join(":".join(row[index].replace(";", "%3B") for row in rows[start:start + 2]) for start in (0, 2))
                        self.assertEqual(actual, value, key)

    def test_tsv_hgvsc_override_is_transcript_specific(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        options.args["outputformat"] = "TSV"
        record = core.Record("chr17\t7675155\t.\tG\tA\t.\tPASS\t.\n", options, None, self.reference)
        variant = record.variants[0]
        variant.addFlag("TRANSCRIPT", "NM_ONE:NM_TWO")
        variant.addFlag("GENE", "ONE:TWO")
        variant.addFlag("CSN", "c.1A>G:c.2C>T")
        record.haplotype_hgvsc_override = "NC_000017.11(NM_ONE):c.[1A>G;2C>T]:NC_000017.11(NM_TWO):c.[3A>G;4C>T]"
        header = io.StringIO()
        core.writeHeader(options, "", header, False, "2.0.15")
        columns = header.getvalue().strip().split("\t")
        out = io.StringIO()
        record.output("TSV", out, options, set(), set(), set(), False)
        rows = [dict(zip(columns, line.split("\t"))) for line in out.getvalue().splitlines()]
        self.assertEqual(rows[0]["HGVSC"], "NC_000017.11(NM_ONE):c.[1A>G;2C>T]")
        self.assertEqual(rows[1]["HGVSC"], "NC_000017.11(NM_TWO):c.[3A>G;4C>T]")
        self.assertEqual(rows[0]["CSNALT"], ".")
        self.assertEqual(rows[1]["CSNALT"], ".")

    def test_tsv_decodes_semicolons_without_changing_vcf_or_other_escapes(self):
        options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
        options.args.update(outputformat="TSV", givealt=True, givealtflag=True)
        options.transcript2protein = {"NM_TEST": "NP_TEST"}
        record = core.Record(
            "chr17\t7675155\t.\tG\tA\t.\tPASS\t.\n",
            options, None, self.reference,
        )
        annotations = {
            "TRANSCRIPT": "NM_TEST",
            "GENE": "TEST",
            "TRINFO": "literal%25and%09",
            "CSN": "c.[1A>G%3B2C>T]_p.[(Met1Val%3bAla2Val)]",
            "CSNALT": "c.1A>G_p.[(Met1Val;Ala2Ala%3BAla3Val)]",
            "ALTANN": "c.[1A>G%3b2C>T]",
            "HGVSg": "NC_000017.11:g.[1A>G%3B2C>T]",
            "CAVA_ORIGHAPLOTYPE": "atom1%3Batom2",
            "CAVA_HAPLOTYPE": "atom1%3batom2",
        }
        for key, value in annotations.items():
            record.variants[0].addFlag(key, value)
        record.haplotype_hgvsc_override = "NC_000017.11(NM_TEST):c.[1A>G%3B2C>T]"
        header = io.StringIO()
        core.writeHeader(options, "", header, False, "2.0.15")
        columns = header.getvalue().strip().split("\t")
        vcf_before = io.StringIO()
        record.output("VCF", vcf_before, options, set(), set(), set(), False)

        for stdout in (False, True):
            with self.subTest(stdout=stdout):
                out = io.StringIO()
                with redirect_stdout(out):
                    record.output("TSV", out, options, set(), set(), set(), stdout)
                row = out.getvalue().strip().split("\t")
                self.assertEqual(len(row), len(columns))
                values = dict(zip(columns, row))
                self.assertNotRegex(out.getvalue(), r"%3[bB]")
                for key in ("CSN", "CSNALT", "ALTANN", "CAVA_ORIGHAPLOTYPE", "CAVA_HAPLOTYPE"):
                    self.assertEqual(values[key], re.sub(r"%3[bB]", ";", annotations[key]))
                self.assertEqual(values["TRINFO"], "literal%25and%09")
                self.assertEqual(values["HGVSG"], "NC_000017.11:g.[1A>G;2C>T]")
                self.assertEqual(values["HGVSC"], "NC_000017.11(NM_TEST):c.[1A>G;2C>T]")
                self.assertEqual(values["HGVSP"], "NP_TEST:p.([(Met1Val;Ala2Val)])")

        vcf_after = io.StringIO()
        record.output("VCF", vcf_after, options, set(), set(), set(), False)
        self.assertEqual(vcf_after.getvalue(), vcf_before.getvalue())
        self.assertIn("%3B", vcf_after.getvalue())
        for key, value in annotations.items():
            self.assertEqual(record.variants[0].getFlag(key), value)

    def test_split_adjacent_protein_has_balanced_delimiters(self):
        annotation = haplotype.maybe_build_splitnearby_csn(
            "c.477_486delinsACGT_p.Met159_Thr162delinsMetAlaValAla",
            "159-162", "MAVT", "MAVA",
        )
        self.assertEqual(
            annotation,
            "c.477_486delinsACGT_p.[(Met159Met;Ala160Ala;Val161Val;Thr162Ala)]",
        )
        self.assertEqual(annotation.count("("), annotation.count(")"))
        self.assertEqual(annotation.count("["), annotation.count("]"))

    def test_split_modes_emit_additional_records(self):
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

            with open(inp, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                # Expected canonical protein delins suitable for split-nearby optional output.
                f.write(
                    "chr17\t7675155\tchr17_7675155_G_A;chr17_7675157_G_C\tGCG\tACC\t.\tPASS\t.\n"
                )

            options = core.Options(cfg)
            options.args["parseHaplotype"] = True
            options.args["splitBasedOnProtein"] = True
            options.args["splitadjacentprotein"] = True

            # Write VCF header first, then execute one in-process job.
            outvcf = outprefix + ".vcf"
            with open(outvcf, "w", encoding="utf-8") as out:
                core.writeHeader(
                    options, "\n".join(main.readHeader(inp)), out, False, "2.0.15"
                )

            impactdir = {}
            for i, valuev in enumerate(options.args["impactdef"].split("|")):
                classv = valuev.split(",")
                for c in classv:
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
                1, options, copts, 3, "", set(), set(), set(), impactdir, 1
            )
            job.run()

            with open(outvcf, "r", encoding="utf-8") as f:
                body_lines = [
                    x.strip() for x in f if not x.startswith("#") and x.strip()
                ]

            # canonical plus at least one additional split/supplemental output
            self.assertGreaterEqual(len(body_lines), 2)
            saw_subset = False
            csn_values = set()
            for line in body_lines:
                self.assertIn("CAVA_ORIGHAPLOTYPE=", line)
                self.assertIn("CAVA_HAPLOTYPE=", line)
                info = line.split("\t")[7]
                self.assertTrue(all("=" in field for field in info.split(";")))
                pairs = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
                for key in ("CSN", "CSNALT", "HGVSp"):
                    value = pairs.get(key, "")
                    self.assertEqual(value.count("("), value.count(")"), value)
                    self.assertEqual(value.count("["), value.count("]"), value)
                if pairs.get("CAVA_ORIGHAPLOTYPE") != pairs.get("CAVA_HAPLOTYPE"):
                    saw_subset = True
                if "CSN" in pairs:
                    csn_values.add(pairs["CSN"])
            self.assertTrue(saw_subset or len(csn_values) > 1)

    def test_record_output_contains_haplotype_tags(self):
        with tempfile.TemporaryDirectory() as td:
            line = "chr17\t7675155\tchr17_7675155_G_A;chr17_7675157_G_C\tGCG\tACC\t.\tPASS\t.\n"
            options = core.Options(os.path.join(self.repo_root, "config_template.txt"))
            options.args["inputformat"] = "VCF"
            options.args["outputformat"] = "VCF"
            options.args["prefix"] = False
            options.args["reference"] = os.path.join(
                self.repo_root, "cava", "data", "tmp.GRCh38.fa"
            )

            record = core.Record(line, options, None, self.reference)
            parsed = haplotype.parse_haplotype_row(
                record.chrom,
                record.pos,
                record.ref,
                record.alts[0],
                record.id,
                self.reference,
                1,
                "smoke",
            )
            original_ids = ";".join([a.token for a in parsed.atomic])
            haplotype.add_haplotype_flags(record, original_ids, original_ids)

            outvcf = os.path.join(td, "out.vcf")
            with open(outvcf, "w", encoding="utf-8") as f:
                record.output("VCF", f, options, set(), set(), set(), False)

            with open(outvcf, "r", encoding="utf-8") as f:
                body = f.read()
            self.assertIn("CAVA_ORIGHAPLOTYPE=", body)
            self.assertIn("CAVA_HAPLOTYPE=", body)

    def test_tsv_output_has_haplotype_columns(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = os.path.join(td, "cfg.txt")
            inp = os.path.join(td, "in.vcf")
            outprefix = os.path.join(td, "out")

            with open(cfg, "w", encoding="utf-8") as f:
                f.write("@inputformat=VCF\n")
                f.write("@outputformat=TSV\n")
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

            with open(inp, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                f.write(
                    "chr17\t7675155\tchr17_7675155_G_A;chr17_7675157_G_C\tGCG\tACC\t.\tPASS\t.\n"
                )

            options = core.Options(cfg)
            options.args["parseHaplotype"] = True

            with open(outprefix + ".txt", "w", encoding="utf-8") as out:
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
                1, options, copts, 3, "", set(), set(), set(), impactdir, 1
            )
            job.run()

            with open(outprefix + ".txt", "r", encoding="utf-8") as f:
                lines = [x.strip() for x in f if x.strip()]

            self.assertTrue(
                lines[0].startswith("ID\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tTYPE")
            )
            self.assertIn("CAVA_ORIGHAPLOTYPE", lines[0])
            self.assertIn("CAVA_HAPLOTYPE", lines[0])
            columns = lines[0].split("\t")
            self.assertIn("CSNALT", columns)
            self.assertEqual(len(columns), len(set(columns)))
            for line in lines[1:]:
                self.assertEqual(len(line.split("\t")), len(columns), line)

            body = "\n".join(lines[1:])
            self.assertIn("chr17_7675155_G_A;chr17_7675157_G_C", body)

    def test_parse_haplotype_hgvsc_uses_delins_for_adjacent_and_cis_for_separated(self):
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

            with open(inp, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                # Adjacent atoms in start codon context.
                f.write(
                    "13\t32316461\t13_32316461_A_C;13_32316462_T_C\tAT\tCC\t.\tPASS\t.\n"
                )
                # Separated atoms (1 bp gap) in start codon context.
                f.write(
                    "13\t32316461\t13_32316461_A_C;13_32316463_G_C\tATG\tCTC\t.\tPASS\t.\n"
                )

            options = core.Options(cfg)
            options.args["parseHaplotype"] = True
            options.args["splitBasedOnProtein"] = False
            options.args["splitadjacentprotein"] = False

            outvcf = outprefix + ".vcf"
            with open(outvcf, "w", encoding="utf-8") as out:
                core.writeHeader(
                    options, "\n".join(main.readHeader(inp)), out, False, "2.0.15"
                )

            impactdir = {}
            for i, valuev in enumerate(options.args["impactdef"].split("|")):
                classv = valuev.split(",")
                for c in classv:
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
                1, options, copts, 3, "", set(), set(), set(), impactdir, 2
            )
            job.run()

            seen = {}
            with open(outvcf, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    if ";" not in cols[2]:
                        continue
                    info = cols[7]
                    pairs = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
                    seen[cols[2]] = {
                        "hgvsc": pairs.get("HGVSc", ""),
                        "hgvsg": pairs.get("HGVSg", ""),
                        "csn": pairs.get("CSN", ""),
                    }

            adj_id = "13_32316461_A_C;13_32316462_T_C"
            sep_id = "13_32316461_A_C;13_32316463_G_C"
            self.assertIn(adj_id, seen)
            self.assertIn(sep_id, seen)
            self.assertNotIn("c.[", seen[adj_id]["hgvsc"])
            self.assertIn("delins", seen[adj_id]["hgvsc"])
            self.assertNotIn("g.[", seen[adj_id]["hgvsg"])
            self.assertIn("delins", seen[adj_id]["hgvsg"])

            self.assertIn("c.[", seen[sep_id]["hgvsc"])
            self.assertIn("%3B", seen[sep_id]["hgvsc"])
            self.assertIn("g.[", seen[sep_id]["hgvsg"])
            self.assertIn("%3B", seen[sep_id]["hgvsg"])
            self.assertNotIn("c.[", seen[sep_id]["csn"])


if __name__ == "__main__":
    unittest.main()
