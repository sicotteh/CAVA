import csv
import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / "scripts/catalogs/prepare_reference_fastas.py"
MAPPER = ROOT / "scripts/catalogs/build_mane_grch37_candidates.py"


def write_gtf(path, transcript_id, symbol="GSTT1"):
    text = (
        'chr1\tfixture\ttranscript\t1\t8\t.\t+\t.\tgene_id "G1"; gene_name "{}"; transcript_id "{}";\n'
        'chr1\tfixture\texon\t1\t8\t.\t+\t.\tgene_id "G1"; gene_name "{}"; transcript_id "{}";\n'
        'chr1\tfixture\tCDS\t1\t3\t.\t+\t0\tgene_id "G1"; gene_name "{}"; transcript_id "{}";\n'
        'chr1\tfixture\tCDS\t6\t8\t.\t+\t0\tgene_id "G1"; gene_name "{}"; transcript_id "{}";\n'
    ).format(symbol, transcript_id, symbol, transcript_id, symbol, transcript_id, symbol, transcript_id)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


class ManeGrch37MappingTests(unittest.TestCase):
    def test_refseq_earlier_version_and_ensembl_require_matching_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp)
            (w / "fasta").mkdir()
            for name in ("grch38_refseq.fa", "grch38_ensembl.fa", "grch37.fa"):
                (w / "fasta" / name).write_text(">chr1\nACGTACGTACGT\n", encoding="ascii")
            subprocess.check_call(["python3", str(PREP), "--workspace", str(w)])

            gr38_refseq = w / "gr38_refseq.gtf.gz"
            gr37_refseq = w / "gr37_refseq.gtf.gz"
            gr38_ensembl = w / "gr38_ensembl.gtf.gz"
            gr37_ensembl = w / "gr37_ensembl.gtf.gz"
            write_gtf(gr38_refseq, "NM_000001.2")
            write_gtf(gr37_refseq, "NM_000001.1")
            write_gtf(gr38_ensembl, "ENST_NEW.1")
            write_gtf(gr37_ensembl, "ENST_OLD.1")

            summary = w / "summary.tsv.gz"
            with gzip.open(summary, "wt", encoding="utf-8") as handle:
                handle.write("symbol\tRefSeq_nuc\tEnsembl_nuc\n")
                handle.write("GSTT1\tNM_000001.2\tENST_NEW.1\n")

            refseq_out = w / "refseq.txt"
            ensembl_out = w / "ensembl.txt"
            report = w / "report.tsv"
            subprocess.check_call([
                "python3", str(MAPPER),
                "--mane-summary", str(summary),
                "--mane-refseq-grch38-gtf", str(gr38_refseq),
                "--mane-ensembl-grch38-gtf", str(gr38_ensembl),
                "--ensembl75-grch37-gtf", str(gr37_ensembl),
                "--grch38-refseq-fasta", str(w / "prepared-fasta/grch38_refseq.fa"),
                "--grch38-fasta", str(w / "prepared-fasta/grch38_ensembl.fa"),
                "--grch37-fasta", str(w / "prepared-fasta/grch37.fa"),
                "--grch37-refseq-gtf", str(gr37_refseq),
                "--refseq-output", str(refseq_out),
                "--ensembl-output", str(ensembl_out),
                "--report", str(report),
            ])
            self.assertEqual(refseq_out.read_text().strip(), "NM_000001.1")
            self.assertEqual(ensembl_out.read_text().strip(), "ENST_OLD.1")
            with report.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["refseq_cds_length_signature_match"], "true")
            self.assertEqual(row["refseq_exact_cds_sequence_match"], "true")
            self.assertIn("cds_lengths_and_sequence_match", row["refseq_method"])
            self.assertEqual(row["ensembl_grch37_candidate"], "ENST_OLD.1")


if __name__ == "__main__":
    unittest.main()
