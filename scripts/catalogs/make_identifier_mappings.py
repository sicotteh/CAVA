#!/usr/bin/env python3
"""Create Ensembl-to-RefSeq mapping files for selenocysteine annotation."""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

ENST = re.compile(r"ENST[0-9]+(?:\.[0-9]+)?")
REFSEQ = re.compile(r"NM_[0-9]+(?:\.[0-9]+)?")


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def from_mane_summary(path: Path) -> set[tuple[str, str]]:
    with open_text(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if "ensembl" in line.lower() and "refseq" in line.lower()
        ),
        None,
    )
    if header_index is None:
        raise SystemExit("Could not identify the MANE summary header")
    header_line = lines[header_index].lstrip("#")
    delimiter = "\t" if "\t" in header_line else ","
    headers = header_line.split(delimiter)
    rows = csv.DictReader(
        (line for line in lines[header_index + 1 :] if not line.startswith("#")),
        fieldnames=headers,
        delimiter=delimiter,
    )
    pairs: set[tuple[str, str]] = set()
    for row in rows:
        text = "\t".join(value or "" for value in row.values())
        ensembl = ENST.findall(text)
        refseq = REFSEQ.findall(text)
        for left in ensembl:
            for right in refseq:
                pairs.add((left, right))
    return pairs


def from_gencode_metadata(path: Path) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            ensembl = ENST.findall(line)
            refseq = REFSEQ.findall(line)
            for left in ensembl:
                for right in refseq:
                    pairs.add((left, right))
    return pairs


def filter_pairs(pairs: set[tuple[str, str]], transcript_file: Path | None) -> set[tuple[str, str]]:
    if transcript_file is None:
        return pairs
    selected = {
        line.strip()
        for line in transcript_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }
    selected_base = {value.split(".", 1)[0] for value in selected}
    return {
        pair
        for pair in pairs
        if pair[0] in selected or pair[0].split(".", 1)[0] in selected_base
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("mane", "gencode"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pairs = from_mane_summary(args.input) if args.mode == "mane" else from_gencode_metadata(args.input)
    pairs = filter_pairs(pairs, args.transcripts)
    if not pairs:
        raise SystemExit("No Ensembl-to-RefSeq mappings were found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "".join("{}\t{}\n".format(left, right) for left, right in sorted(pairs)),
        encoding="utf-8",
    )
    print("{}: {} mappings".format(args.output, len(pairs)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
