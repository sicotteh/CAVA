import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import gzip

ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / "scripts" / "catalogs" / "prepare_reference_fastas.py"
spec = importlib.util.spec_from_file_location("prepare_reference_fastas", PREP)
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)


class CatalogSequenceValidationTests(unittest.TestCase):
    def test_plus_and_minus_catalog_coordinates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nAAAACCCCGGGGTTTT\n", encoding="ascii")
            prep.build_fai(fasta)
            rna = root / "rna.fa"
            rna.write_text(">TXP\nAAAAGGGG\n>TXM\nAAAAGGGG\n", encoding="ascii")
            catalog = root / "catalog.gz"
            plus = "TXP\tGENEP\tGENEP\t+/12bp/2/8bp/1\tchr1\t1\t0\t12\t0\t1\t12\t0\t4\t8\t12\n"
            minus = "TXM\tGENEM\tGENEM\t-/12bp/2/8bp/1\tchr1\t-1\t4\t16\t0\t16\t5\t12\t16\t4\t8\n"
            with gzip.open(catalog, "wt", encoding="utf-8") as handle:
                handle.write(plus)
                handle.write(minus)
            report = root / "report.tsv"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "catalogs" / "validate_catalog_sequences.py"),
                    "--catalog",
                    str(catalog),
                    "--genome-fasta",
                    str(fasta),
                    "--rna-fasta",
                    str(rna),
                    "--report",
                    str(report),
                ],
                cwd=ROOT / "scripts" / "catalogs",
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            text = report.read_text(encoding="utf-8")
            self.assertIn("TXP", text)
            self.assertIn("TXM", text)
            self.assertNotIn("sequence_mismatch", text)


if __name__ == "__main__":
    unittest.main()
