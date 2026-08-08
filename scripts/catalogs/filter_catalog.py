#!/usr/bin/env python3
"""Filter a CAVA catalog/map to validated transcript IDs and reindex it."""
from __future__ import annotations

import argparse
import gzip
import os
import tempfile
from pathlib import Path

import pysam


def read_ids(path: Path) -> set[str]:
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--transcript-map", type=Path, required=True)
    parser.add_argument("--allowed", type=Path, required=True)
    parser.add_argument("--output-catalog", type=Path, required=True)
    parser.add_argument("--output-map", type=Path, required=True)
    parser.add_argument("--require-gene", action="append", default=[])
    args = parser.parse_args()

    allowed = read_ids(args.allowed)
    if not allowed:
        raise SystemExit("The allowed transcript list is empty")
    args.output_catalog.parent.mkdir(parents=True, exist_ok=True)
    args.output_map.parent.mkdir(parents=True, exist_ok=True)

    found_ids: set[str] = set()
    found_genes: set[str] = set()
    with tempfile.TemporaryDirectory(prefix="cava-filter-") as temporary:
        plain = Path(temporary) / "catalog.tsv"
        opener = gzip.open if args.catalog.suffix == ".gz" else open
        with opener(args.catalog, "rt", encoding="utf-8", errors="replace") as source, plain.open(
            "w", encoding="utf-8", newline="\n"
        ) as target:
            for line in source:
                if not line or line.startswith("#"):
                    target.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 13:
                    raise SystemExit("Malformed CAVA catalog line: {}".format(line[:200]))
                if fields[0] not in allowed:
                    continue
                target.write(line)
                found_ids.add(fields[0])
                found_genes.add(fields[1].upper())
        pysam.tabix_compress(str(plain), str(args.output_catalog), force=True)
        pysam.tabix_index(
            str(args.output_catalog),
            seq_col=4,
            start_col=6,
            end_col=7,
            zerobased=True,
            meta_char="#",
            force=True,
        )

    with args.transcript_map.open("r", encoding="utf-8", errors="replace") as source, args.output_map.open(
        "w", encoding="utf-8", newline="\n"
    ) as target:
        for line in source:
            if not line or line.startswith("#") or line.startswith("GENEID"):
                target.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] in found_ids:
                target.write(line)

    missing_genes = sorted({gene.upper() for gene in args.require_gene} - found_genes)
    if missing_genes:
        raise SystemExit(
            "Required genes are absent after filtering: {}".format(", ".join(missing_genes))
        )
    missing_ids = allowed - found_ids
    print(
        "{}: retained={} requested_but_not_built={}".format(
            args.output_catalog, len(found_ids), len(missing_ids)
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
