#!/usr/bin/env python3
"""Create deterministic transcript lists from GTF attributes."""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

TRANSCRIPT_RE = re.compile(r'(?:^|;\s*)transcript_id\s+"([^"]+)"')
TAG_MANE_RE = re.compile(r'(?:^|;\s*)tag\s+"MANE (?:Select|Plus Clinical)"')


def ids_from_gtf(path: Path, *, mane_only: bool = False) -> list[str]:
    identifiers: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            attributes = fields[8]
            if mane_only and not TAG_MANE_RE.search(attributes):
                continue
            match = TRANSCRIPT_RE.search(attributes)
            if match:
                identifiers.add(match.group(1))
    return sorted(identifiers)


def write_ids(path: Path, identifiers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(identifier + "\n" for identifier in identifiers), encoding="utf-8")
    print(f"{path}: {len(identifiers)} transcripts")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, default=Path("catalog-build"))
    args = parser.parse_args()
    w = args.workspace
    jobs = [
        (w / "gtf/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz", w / "transcripts/MANE_1.5_GRCh38_REFSEQ.txt", False),
        (w / "gtf/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", w / "transcripts/MANE_1.5_GRCh38_ENSEMBL.txt", False),
        (w / "gtf/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz", w / "transcripts/Ensembl_116_GRCh38.txt", False),
        (w / "gtf/Homo_sapiens.GRCh37.75.gtf.gz", w / "transcripts/Ensembl_75_GRCh37.txt", False),
        (w / "gtf/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz", w / "transcripts/GENCODE_50_GRCh38.txt", False),
        (w / "gtf/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz", w / "transcripts/MANE_1.5_GRCh37_REFSEQ_candidates.txt", True),
    ]
    for source, target, mane_only in jobs:
        if not source.is_file():
            raise SystemExit(f"Missing source GTF: {source}")
        write_ids(target, ids_from_gtf(source, mane_only=mane_only))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
