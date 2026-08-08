import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / "scripts/catalogs/prepare_reference_fastas.py"
VALIDATE = ROOT / "scripts/catalogs/validate_gtf_sequences.py"
NORMALIZE = ROOT / "scripts/catalogs/normalize_gtf_contigs.py"


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class SequenceValidationTests(unittest.TestCase):
    def test_one_based_plus_and_minus_transcripts(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp)
            (w / "fasta").mkdir()
            sequence = ("ACGT" * 30)[:100]
            (w / "fasta/test.fa").write_text(">chr1\n" + sequence + "\n", encoding="ascii")
            subprocess.check_call(["python3", str(PREP), "--workspace", str(w)])
            gtf = w / "test.gtf"
            gtf.write_text(
                "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_name \"PLUS\"; transcript_id \"TXP.1\";\n"
                "chr1\ttest\texon\t9\t12\t.\t+\t.\tgene_name \"PLUS\"; transcript_id \"TXP.1\";\n"
                "chr1\ttest\texon\t21\t24\t.\t-\t.\tgene_name \"MINUS\"; transcript_id \"TXM.1\";\n"
                "chr1\ttest\texon\t29\t32\t.\t-\t.\tgene_name \"MINUS\"; transcript_id \"TXM.1\";\n",
                encoding="utf-8",
            )
            plus = sequence[0:4] + sequence[8:12]
            minus = rc(sequence[28:32]) + rc(sequence[20:24])
            rna = w / "rna.fa"
            rna.write_text(">TXP.1\n{}\n>TXM.1\n{}\n".format(plus, minus), encoding="ascii")
            tx = w / "tx.txt"; tx.write_text("TXP.1\nTXM.1\n", encoding="ascii")
            report = w / "report.tsv"; accepted = w / "accepted.txt"
            subprocess.check_call([
                "python3", str(VALIDATE), "--gtf", str(gtf),
                "--genome-fasta", str(w / "prepared-fasta/test.fa"),
                "--rna-fasta", str(rna), "--transcripts", str(tx),
                "--report", str(report), "--accepted-transcripts", str(accepted),
            ])
            with report.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["status"] for row in rows}, {"ok"})
            self.assertEqual(set(accepted.read_text().split()), {"TXP.1", "TXM.1"})

    def test_assembly_report_alias_resolves_alt_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp)
            (w / "fasta").mkdir()
            (w / "fasta/reference.fa").write_text(">NC_TEST.1\nACGTACGT\n", encoding="ascii")
            (w / "fasta/test_assembly_report.txt").write_text(
                "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
                "ALT_TEST\talt-scaffold\tna\tna\tGB_TEST.1\t=\tNC_TEST.1\tPrimary Assembly\t8\tchrALT_TEST\n",
                encoding="utf-8",
            )
            subprocess.check_call(["python3", str(PREP), "--workspace", str(w)])
            gtf = w / "test.gtf"
            gtf.write_text(
                "chrALT_TEST\ttest\texon\t1\t8\t.\t+\t.\tgene_name \"GSTT1\"; transcript_id \"TXALT.1\";\n",
                encoding="utf-8",
            )
            rna = w / "rna.fa"; rna.write_text(">TXALT.1\nACGTACGT\n", encoding="ascii")
            tx = w / "tx.txt"; tx.write_text("TXALT.1\n", encoding="ascii")
            report = w / "report.tsv"; accepted = w / "accepted.txt"
            subprocess.check_call([
                "python3", str(VALIDATE), "--gtf", str(gtf),
                "--genome-fasta", str(w / "prepared-fasta/reference.fa"),
                "--rna-fasta", str(rna), "--transcripts", str(tx),
                "--report", str(report), "--accepted-transcripts", str(accepted),
                "--require-gene", "GSTT1",
            ])
            self.assertEqual(accepted.read_text().strip(), "TXALT.1")

    def test_normalized_gtf_uses_exact_fasta_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp)
            (w / "fasta").mkdir()
            (w / "fasta/reference.fa").write_text(">NC_TEST.1\nACGTACGT\n", encoding="ascii")
            (w / "fasta/test_assembly_report.txt").write_text(
                "ALT_TEST\talt-scaffold\tna\tna\tGB_TEST.1\t=\tNC_TEST.1\tPrimary Assembly\t8\tchrALT_TEST\n",
                encoding="utf-8",
            )
            subprocess.check_call(["python3", str(PREP), "--workspace", str(w)])
            gtf = w / "input.gtf"
            gtf.write_text(
                "chrALT_TEST\ttest\tgene\t1\t8\t.\t+\t.\tgene_id \"G1\"; gene_name \"GSTT1\";\n"
                "chrALT_TEST\ttest\ttranscript\t1\t8\t.\t+\t.\tgene_id \"G1\"; gene_name \"GSTT1\"; transcript_id \"TXALT.1\";\n"
                "chrALT_TEST\ttest\texon\t1\t8\t.\t+\t.\tgene_id \"G1\"; gene_name \"GSTT1\"; transcript_id \"TXALT.1\";\n",
                encoding="utf-8",
            )
            selected = w / "selected.txt"; selected.write_text("TXALT.1\n", encoding="utf-8")
            output = w / "normalized.gtf.gz"; report = w / "mapping.tsv"
            subprocess.check_call([
                "python3", str(NORMALIZE), "--gtf", str(gtf),
                "--genome-fasta", str(w / "prepared-fasta/reference.fa"),
                "--transcripts", str(selected), "--output", str(output),
                "--mapping-report", str(report),
            ])
            import gzip
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                lines = [line for line in handle if not line.startswith("#")]
            self.assertTrue(lines)
            self.assertTrue(all(line.startswith("NC_TEST.1\t") for line in lines))
            self.assertIn("chrALT_TEST\tNC_TEST.1", report.read_text(encoding="utf-8"))

    def test_duplicate_accession_selects_exact_alt_locus(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp)
            (w / "fasta").mkdir()
            # Primary locus is an exact duplicate in structure but has the wrong
            # sequence. The alternate locus is the authoritative RNA match.
            (w / "fasta/reference.fa").write_text(
                ">chr22\nAAAAAAAA\n>chr22_KI270879v1_alt\nACGTACGT\n",
                encoding="ascii",
            )
            subprocess.check_call(["python3", str(PREP), "--workspace", str(w)])
            gtf = w / "gstt1.gtf"
            attrs = 'gene_name "GSTT1"; transcript_id "NM_000853.4";'
            gtf.write_text(
                "chr22\ttest\texon\t1\t8\t.\t+\t.\t" + attrs + "\n"
                "chr22_KI270879v1_alt\ttest\texon\t1\t8\t.\t+\t.\t" + attrs + "\n",
                encoding="utf-8",
            )
            rna = w / "rna.fa"
            rna.write_text(">NM_000853.4\nACGTACGT\n", encoding="ascii")
            tx = w / "tx.txt"
            tx.write_text("NM_000853.4\n", encoding="ascii")
            report = w / "report.tsv"
            accepted = w / "accepted.txt"
            filtered = w / "filtered.gtf.gz"
            subprocess.check_call([
                "python3", str(VALIDATE), "--gtf", str(gtf),
                "--genome-fasta", str(w / "prepared-fasta/reference.fa"),
                "--rna-fasta", str(rna), "--transcripts", str(tx),
                "--report", str(report), "--accepted-transcripts", str(accepted),
                "--filtered-gtf", str(filtered), "--require-gene", "GSTT1",
            ])
            import gzip
            with gzip.open(filtered, "rt", encoding="utf-8") as handle:
                selected_lines = [line for line in handle if not line.startswith("#")]
            self.assertEqual(len(selected_lines), 1)
            self.assertTrue(selected_lines[0].startswith("chr22_KI270879v1_alt\t"))
            with report.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            status_by_contig = {row["contig"]: row["status"] for row in rows}
            self.assertEqual(status_by_contig["chr22"], "sequence_mismatch")
            self.assertEqual(status_by_contig["chr22_KI270879v1_alt"], "ok")


if __name__ == "__main__":
    unittest.main()
