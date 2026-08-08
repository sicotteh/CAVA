#!/usr/bin/env python3
"""Validate cava_catalogs.tsv and optionally materialized Git LFS payloads."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path, PurePosixPath

FILE_GROUPS = (
    ("lfs_path", "file_name", "sha256", "lfs_oid_sha256", "lfs_size"),
    ("index_path", "index_file_name", "index_sha256", "index_lfs_oid_sha256", "index_lfs_size"),
    (
        "transcript_map_path",
        "transcript_map_file_name",
        "transcript_map_sha256",
        "transcript_map_lfs_oid_sha256",
        "transcript_map_lfs_size",
    ),
    (
        "selenocysteine_path",
        "selenocysteine_file_name",
        "selenocysteine_sha256",
        "selenocysteine_lfs_oid_sha256",
        "selenocysteine_lfs_size",
    ),
)
REQUIRED = [
    "catalog_name",
    "file_name",
    "Build",
    "source",
    "source_version",
    "transcript_ID_type",
    "repository",
    "ref",
    "status",
]
for group in FILE_GROUPS:
    REQUIRED.extend(group)
HEX = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_path(value: str) -> bool:
    path = PurePosixPath(value)
    return bool(value) and not path.is_absolute() and ".." not in path.parts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("cava_catalogs.tsv"))
    parser.add_argument("--check-files", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    manifest = args.manifest if args.manifest.is_absolute() else repo / args.manifest
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    errors: list[str] = []
    for column in REQUIRED:
        if column not in fieldnames:
            errors.append("missing column: {}".format(column))
    names: set[str] = set()
    for number, row in enumerate(rows, 2):
        name = row.get("catalog_name", "")
        prefix = "row {} ({})".format(number, name or "unnamed")
        if not name or name in names:
            errors.append("{}: catalog_name is empty or duplicated".format(prefix))
        names.add(name)
        if row.get("Build") not in {"GRCh38", "GRCh37"}:
            errors.append("{}: Build must be GRCh38 or GRCh37".format(prefix))
        if row.get("transcript_ID_type") not in {"ENSEMBL", "REFSEQ"}:
            errors.append(
                "{}: transcript_ID_type must be ENSEMBL or REFSEQ".format(prefix)
            )
        if row.get("source", "").upper() != "MANE" and row.get("transcript_ID_type") == "REFSEQ":
            errors.append("{}: only MANE rows may use REFSEQ IDs".format(prefix))
        status = row.get("status", "")
        if status not in {"planned", "experimental", "published", "retired"}:
            errors.append("{}: unsupported status {!r}".format(prefix, status))

        for path_key, name_key, hash_key, oid_key, size_key in FILE_GROUPS:
            relative = row.get(path_key, "")
            if not safe_path(relative):
                errors.append("{}: unsafe {}: {!r}".format(prefix, path_key, relative))
                continue
            if Path(relative).name != row.get(name_key):
                errors.append("{}: {} does not match {}".format(prefix, name_key, path_key))
            materialized = status == "published" or (
                status == "experimental"
                and HEX.fullmatch(row.get(hash_key, "")) is not None
                and row.get(size_key, "").isdigit()
            )
            if status == "published":
                if not HEX.fullmatch(row.get(hash_key, "")):
                    errors.append("{}: invalid {}".format(prefix, hash_key))
                if not HEX.fullmatch(row.get(oid_key, "")):
                    errors.append("{}: invalid {}".format(prefix, oid_key))
                if row.get(hash_key) != row.get(oid_key):
                    errors.append("{}: {} and {} differ".format(prefix, hash_key, oid_key))
                if not row.get(size_key, "").isdigit():
                    errors.append("{}: invalid {}".format(prefix, size_key))
            if args.check_files and materialized:
                path = repo / relative
                if not path.is_file():
                    errors.append("{}: missing worktree file {}".format(prefix, path))
                    continue
                expected_hash = row.get(hash_key, "")
                expected_size = row.get(size_key, "")
                if HEX.fullmatch(expected_hash) and sha256(path) != expected_hash:
                    errors.append("{}: checksum mismatch for {}".format(prefix, path))
                if expected_size.isdigit() and path.stat().st_size != int(expected_size):
                    errors.append("{}: size mismatch for {}".format(prefix, path))
    if not rows:
        errors.append("manifest has no rows")
    if errors:
        for error in errors:
            print("ERROR:", error)
        return 1
    print("{}: {} catalog rows are valid".format(manifest, len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
