#!/usr/bin/env python3
"""Copy one built catalog and its sidecars into LFS paths and update its row."""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

FILE_GROUPS = (
    ("catalog", "lfs_path", "file_name", "sha256", "lfs_oid_sha256", "lfs_size"),
    ("index", "index_path", "index_file_name", "index_sha256", "index_lfs_oid_sha256", "index_lfs_size"),
    (
        "transcript_map",
        "transcript_map_path",
        "transcript_map_file_name",
        "transcript_map_sha256",
        "transcript_map_lfs_oid_sha256",
        "transcript_map_lfs_size",
    ),
    (
        "selenocysteine",
        "selenocysteine_path",
        "selenocysteine_file_name",
        "selenocysteine_sha256",
        "selenocysteine_lfs_oid_sha256",
        "selenocysteine_lfs_size",
    ),
)
REQUIRED = {"catalog_name", "status"}
for _, path_key, name_key, hash_key, oid_key, size_key in FILE_GROUPS:
    REQUIRED.update((path_key, name_key, hash_key, oid_key, size_key))


def digest(path: Path) -> tuple[str, int]:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest(), path.stat().st_size


def copy_atomic(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == target.resolve():
        return
    temporary = target.with_name(target.name + ".part")
    shutil.copy2(source, temporary)
    os.replace(temporary, target)


def check_lfs_attribute(repo: Path, relative: str) -> None:
    result = subprocess.run(
        ["git", "check-attr", "filter", "--", relative],
        cwd=repo,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    if not result.stdout.rstrip().endswith(": lfs"):
        raise SystemExit(
            "{} is not tracked by Git LFS; run setup_git_lfs.sh first".format(relative)
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog_name")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("cava_catalogs.tsv"))
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--transcript-map", type=Path, required=True)
    parser.add_argument("--selenocysteine", type=Path, required=True)
    parser.add_argument(
        "--status", choices=("published", "experimental"), default="published"
    )
    parser.add_argument("--no-copy", action="store_true")
    args = parser.parse_args()

    repo = args.repo.resolve()
    manifest = args.manifest if args.manifest.is_absolute() else repo / args.manifest
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    missing = sorted(REQUIRED - set(fieldnames))
    if missing:
        raise SystemExit("Manifest is missing columns: " + ", ".join(missing))

    matches = [row for row in rows if row["catalog_name"] == args.catalog_name]
    if len(matches) != 1:
        raise SystemExit(
            "Expected exactly one row named {!r}; found {}".format(
                args.catalog_name, len(matches)
            )
        )
    row = matches[0]
    sources = {
        "catalog": args.catalog.resolve(),
        "index": args.index.resolve(),
        "transcript_map": args.transcript_map.resolve(),
        "selenocysteine": args.selenocysteine.resolve(),
    }

    for label, source in sources.items():
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit("Missing or empty --{} file: {}".format(label.replace("_", "-"), source))

    for label, path_key, name_key, hash_key, oid_key, size_key in FILE_GROUPS:
        target = repo / row[path_key]
        check_lfs_attribute(repo, row[path_key])
        if not args.no_copy:
            copy_atomic(sources[label], target)
        if not target.is_file():
            raise SystemExit("Expected LFS worktree file does not exist: {}".format(target))
        sha, size = digest(target)
        row[name_key] = target.name
        row[hash_key] = sha
        # A Git LFS object ID is the SHA-256 of materialized content.
        row[oid_key] = sha
        row[size_key] = str(size)
        print("{}: sha256={} size={}".format(target, sha, size))

    row["status"] = args.status
    fd, temporary_name = tempfile.mkstemp(prefix=manifest.name + ".", dir=manifest.parent)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(
                out,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, manifest)
    finally:
        temporary.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
