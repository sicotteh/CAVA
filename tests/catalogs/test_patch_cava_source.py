import importlib.util
import tempfile
import tomllib
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2] / "scripts" / "catalogs" / "patch_cava_source.py"
spec = importlib.util.spec_from_file_location("patch_cava_source", MODULE)
patcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(patcher)


class PatchCavaSourceTests(unittest.TestCase):
    def make_repo(self, root: Path) -> None:
        (root / "cava" / "ensembldb").mkdir(parents=True)
        (root / "cava" / "utils").mkdir(parents=True)
        (root / "pyproject.toml").write_text(
            "[project]\nname = \"CAVA\"\n\n[tool.setuptools]\ninclude-package-data = true\n",
            encoding="utf-8",
        )
        (root / "cava" / "MANE.py").write_text(
            "from optparse import OptionParser\n"
            "from cava.ensembldb import mane_db_prep as main\n"
            "import os\n"
            "parser = OptionParser()\n"
            "parser.add_option(\"--no_hg19\",  action='store_false', default=True, dest='no_hg19', help=\"Set this to skip hg19 builds\")\n"
            "(options, args) = parser.parse_args()\n\n"
            "options.select = False\n",
            encoding="utf-8",
        )
        (root / "cava" / "ensembldb" / "mane_db_prep.py").write_text(
            "import os\n"
            "failed_conversions = dict()\n"
            "failed_conversions['ENST'] = set()\n\n\n"
            "def getValue(tags, tag):\n"
            "    return None\n\n"
            "def parse_GTF(filename='', genesdata=None):\n"
            "    for line in []:\n"
            "        cols = line.split('\\t')\n\n"
            "        # Only consider transcripts on the following chromosomes\n"
            "        if cols[0] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',\n"
            "                           '18', '19', '20', '21', '22', '23', 'MT', 'X', 'Y']: continue\n\n"
            "        # Consider only certain types of lines\n"
            "        if cols[2] not in ['exon', 'transcript', 'start_codon', 'stop_codon']: continue\n\n"
            "        # Annotation tags\n"
            "        tags = cols[8].split(';')\n"
            "    return None, None, None, genesdata\n\n"
            "def parse_gtf_loop(source_compressed_gtf, options, genesdata, transIDs):\n"
            "    transcript, prevenst, first, genesdata = parse_GTF(filename=source_compressed_gtf,\n"
            "                                                       genesdata=genesdata)\n"
            "    return 0\n",
            encoding="utf-8",
        )
        (root / "cava" / "utils" / "data.py").write_text(
            "import os\nimport re\n"
            "def f(options):\n"
            "        if False:\n"
            "            fid = None\n"
            "        else:\n"
            "            fid = _open_ensembldb_resource(\"SECIS_in_refseq_pos.txt\")\n",
            encoding="utf-8",
        )

    def test_patch_is_idempotent_and_adds_gstt1_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_repo(root)
            self.assertTrue(patcher.patch_pyproject(root / "pyproject.toml"))
            self.assertTrue(patcher.patch_mane_cli(root / "cava" / "MANE.py"))
            self.assertTrue(
                patcher.patch_mane_builder(root / "cava" / "ensembldb" / "mane_db_prep.py")
            )
            self.assertTrue(
                patcher.patch_selenodata_discovery(root / "cava" / "utils" / "data.py")
            )
            self.assertFalse(patcher.patch_pyproject(root / "pyproject.toml"))
            self.assertFalse(patcher.patch_mane_cli(root / "cava" / "MANE.py"))
            self.assertFalse(
                patcher.patch_mane_builder(root / "cava" / "ensembldb" / "mane_db_prep.py")
            )
            self.assertFalse(
                patcher.patch_selenodata_discovery(root / "cava" / "utils" / "data.py")
            )

            mane = (root / "cava" / "MANE.py").read_text(encoding="utf-8")
            self.assertIn("--include-alt-gene", mane)
            self.assertIn("--include-alt-transcript", mane)
            builder = (root / "cava" / "ensembldb" / "mane_db_prep.py").read_text(
                encoding="utf-8"
            )
            self.assertIn("def _mane_record_is_selected", builder)
            self.assertIn("options=options", builder)
            self.assertIn("_mane_transcript_base", builder)
            runtime = (root / "cava" / "utils" / "data.py").read_text(encoding="utf-8")
            self.assertIn("adjacent_selenofile", runtime)
            with (root / "pyproject.toml").open("rb") as handle:
                project = tomllib.load(handle)
            self.assertEqual(
                project["project"]["scripts"]["cava_data"],
                "cava.cava_data:main",
            )
            self.assertEqual(
                project["tool"]["setuptools"]["data-files"]["share/cava"],
                ["config_template.txt", "cava_catalogs.tsv"],
            )

    def test_mane_alt_selector_is_opt_in(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_repo(root)
            path = root / "cava" / "ensembldb" / "mane_db_prep.py"
            patcher.patch_mane_builder(path)
            source = path.read_text(encoding="utf-8")
            start = source.index("PRIMARY_MANE_CONTIGS")
            end = source.index("def getValue", start)
            namespace = {"re": __import__("re")}

            def get_value(tags, tag):
                for item in tags:
                    item = item.strip()
                    if item.startswith(tag + " "):
                        return item.split('"', 2)[1]
                return None

            namespace["getValue"] = get_value
            exec(source[start:end], namespace)
            selected = namespace["_mane_record_is_selected"]
            tags = [
                'gene_id "GSTT1"',
                'gene_name "GSTT1"',
                'transcript_id "NM_000853.4_1"',
            ]

            class Options:
                include_alt_genes = set()
                include_alt_transcripts = set()

            options = Options()
            self.assertTrue(selected("22", tags, options))
            self.assertFalse(selected("chr22_KI270879v1_alt", tags, options))
            options.include_alt_genes = {"GSTT1"}
            self.assertTrue(selected("chr22_KI270879v1_alt", tags, options))
            options.include_alt_genes = set()
            options.include_alt_transcripts = {"NM_000853.4"}
            self.assertTrue(selected("chr22_KI270879v1_alt", tags, options))
            other = [
                'gene_name "OTHER"',
                'transcript_id "NM_999999.1"',
            ]
            self.assertFalse(selected("chr1_alt", other, options))


if __name__ == "__main__":
    unittest.main()
