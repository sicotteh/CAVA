import contextlib
import csv
import hashlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

MODULE_PATH = Path(__file__).resolve().parents[2] / "cava" / "cava_data.py"
spec = importlib.util.spec_from_file_location("cava_data_under_test", str(MODULE_PATH))
cava_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cava_data)


def digest(data):
    return hashlib.sha256(data).hexdigest()


class CavaDataTests(unittest.TestCase):
    def write_manifest(self, root, status="published", build="GRCh37"):
        payloads = {
            "catalog": b"catalog\n",
            "index": b"index\n",
            "transcript_map": b"gene\tgene\ttx\tprotein\n",
            "selenocysteine": b"gene\tSECIS_maxstart\n",
        }
        path = root / "manifest.tsv"
        fields = [
            "catalog_name", "file_name", "Build", "source", "source_version",
            "transcript_ID_type", "repository", "ref", "lfs_path",
            "index_file_name", "index_path", "transcript_map_file_name",
            "transcript_map_path", "selenocysteine_file_name",
            "selenocysteine_path", "sha256", "index_sha256",
            "transcript_map_sha256", "selenocysteine_sha256",
            "lfs_oid_sha256", "lfs_size", "index_lfs_oid_sha256",
            "index_lfs_size", "transcript_map_lfs_oid_sha256",
            "transcript_map_lfs_size", "selenocysteine_lfs_oid_sha256",
            "selenocysteine_lfs_size", "status",
        ]
        row = {
            "catalog_name": "example", "file_name": "example.gz", "Build": build,
            "source": "MANE", "source_version": "1.5", "transcript_ID_type": "REFSEQ",
            "repository": "fake", "ref": "master", "lfs_path": "catalogs/files/example.gz",
            "index_file_name": "example.gz.tbi", "index_path": "catalogs/files/example.gz.tbi",
            "transcript_map_file_name": "example.txt", "transcript_map_path": "catalogs/files/example.txt",
            "selenocysteine_file_name": "example.cesis", "selenocysteine_path": "catalogs/files/example.cesis",
            "status": status,
        }
        for label, prefix in (("catalog", ""), ("index", "index_"), ("transcript_map", "transcript_map_"), ("selenocysteine", "selenocysteine_")):
            value = payloads[label]
            row[prefix + "sha256"] = digest(value)
            row[prefix + "lfs_oid_sha256"] = digest(value)
            row[prefix + "lfs_size"] = str(len(value))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        return path, payloads

    def test_hg19_alias_lists_grch37(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, _ = self.write_manifest(Path(tmp))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = cava_data.main(["--manifest", str(manifest), "list", "--build", "hg19"])
            self.assertEqual(code, 0)
            self.assertIn("example", output.getvalue())

    def test_install_fetches_all_sidecars_and_validates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, payloads = self.write_manifest(root)
            catalog = cava_data.read_manifest(str(manifest))[0]
            destination = root / "installed"

            def fake_clone(selected, checkout):
                for repository_path, output_name, _, _, _ in selected.install_files:
                    target = checkout / repository_path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    label = {
                        "example.gz": "catalog",
                        "example.gz.tbi": "index",
                        "example.txt": "transcript_map",
                        "example.cesis": "selenocysteine",
                    }[output_name]
                    target.write_bytes(payloads[label])

            commands = []
            def fake_run(command, **kwargs):
                commands.append(command)
                return mock.Mock()

            with mock.patch.object(cava_data, "_require_git_lfs"), \
                 mock.patch.object(cava_data, "_clone_pointer_checkout", side_effect=fake_clone), \
                 mock.patch.object(cava_data, "_run", side_effect=fake_run):
                paths = cava_data.install_catalog(catalog, destination)
            self.assertEqual(len(paths), 4)
            for name in ("example.gz", "example.gz.tbi", "example.txt", "example.cesis"):
                self.assertTrue((destination / name).is_file())
            include = "--include=catalogs/files/example.gz,catalogs/files/example.gz.tbi,catalogs/files/example.txt,catalogs/files/example.cesis"
            self.assertTrue(any(include in command for command in commands))

    def test_config_changes_only_ensembl_setting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, payloads = self.write_manifest(root)
            install = root / "data"
            install.mkdir()
            for name, label in (
                ("example.gz", "catalog"),
                ("example.gz.tbi", "index"),
                ("example.txt", "transcript_map"),
                ("example.cesis", "selenocysteine"),
            ):
                (install / name).write_bytes(payloads[label])
            template = root / "config_template.txt"
            template.write_text(
                "@inputformat = VCF\n@ensembl = old.db\n@outputformat = VCF\n",
                encoding="utf-8",
            )
            output = root / "config.txt"
            code = cava_data.main([
                "--manifest", str(manifest), "config", "example", str(install),
                "--template", str(template), "--output", str(output),
            ])
            self.assertEqual(code, 0)
            rendered = output.read_text(encoding="utf-8")
            self.assertIn("@inputformat = VCF", rendered)
            self.assertIn("@outputformat = VCF", rendered)
            self.assertIn("@ensembl = " + str((install / "example.gz").resolve()), rendered)
            self.assertNotIn("old.db", rendered)


if __name__ == "__main__":
    unittest.main()
