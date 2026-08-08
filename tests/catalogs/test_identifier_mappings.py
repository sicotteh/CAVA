import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2] / "scripts" / "catalogs" / "make_identifier_mappings.py"
spec = importlib.util.spec_from_file_location("make_identifier_mappings", MODULE)
mappings = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mappings)


class IdentifierMappingTests(unittest.TestCase):
    def test_mane_summary_header_and_versioned_accessions(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "MANE.summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("# MANE synthetic fixture\n")
                handle.write("#symbol\tRefSeq_nuc\tEnsembl_nuc\tMANE_status\n")
                handle.write("GSTT1\tNM_000853.4\tENST00000612885.4\tMANE Select\n")
            self.assertEqual(
                mappings.from_mane_summary(path),
                {("ENST00000612885.4", "NM_000853.4")},
            )

    def test_gencode_metadata_is_filterable_versionlessly(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.gz"
            with gzip.open(metadata, "wt", encoding="utf-8") as handle:
                handle.write("ENST000001.9\tNM_000001.2\n")
                handle.write("ENST000002.3\tNM_000002.4\n")
            selected = root / "selected.txt"
            selected.write_text("ENST000002\n", encoding="utf-8")
            pairs = mappings.filter_pairs(
                mappings.from_gencode_metadata(metadata), selected
            )
            self.assertEqual(pairs, {("ENST000002.3", "NM_000002.4")})


if __name__ == "__main__":
    unittest.main()
