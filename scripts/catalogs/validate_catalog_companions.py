#!/usr/bin/env python3
"""Validate companion files associated with a final CAVA catalog.

Rules:
1) Every catalog transcript must exist in the companion transcript map.
2) Missing protein IDs are only a QC failure when the transcript has CDS in
    the catalog (protein-coding expectation). For non-coding transcripts (no CDS
    in catalog), missing protein IDs are reported as informational only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def catalog_transcripts(path: Path) -> tuple[set[str], set[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    values: set[str] = set()
    has_cds: set[str] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                raise SystemExit(f"Malformed catalog line in {path}: {line[:200]}")
            transcript = fields[0]
            values.add(transcript)
            try:
                coding_start_rel = int(fields[8])
                coding_start = int(fields[9])
                coding_end = int(fields[10])
            except ValueError:
                coding_start_rel = coding_start = coding_end = -1
            if coding_start_rel > 0 and coding_start > 0 and coding_end > 0 and coding_end >= coding_start:
                has_cds.add(transcript)
    if not values:
        raise SystemExit(f"No transcript rows found in catalog: {path}")
    return values, has_cds


def read_transcript_map(path: Path) -> tuple[set[str], set[str], dict[str, str]]:
    transcripts: set[str] = set()
    proteins_present: set[str] = set()
    genes_by_tx: dict[str, str] = {}
    expected_header = ["GENEID", "SYMBOL", "Transcript", "Protein"]
    saw_header = False
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            first = row[0].strip()
            normalized = [cell.strip() for cell in row]

            if first.startswith("#"):
                candidate = [first.lstrip("#").strip()] + normalized[1:4]
                if len(candidate) >= 4 and candidate[:4] == expected_header:
                    saw_header = True
                continue

            if len(normalized) >= 4 and normalized[:4] == expected_header:
                saw_header = True
                continue

            if len(normalized) < 4:
                raise SystemExit(f"Malformed transcript map row in {path}: {row}")

            transcript = normalized[2]
            protein = normalized[3]
            if not transcript:
                raise SystemExit(f"Blank transcript ID in transcript map row: {row}")
            transcripts.add(transcript)
            genes_by_tx[transcript] = normalized[1] or normalized[0].lstrip("#")
            if protein:
                proteins_present.add(transcript)

    if not transcripts:
        raise SystemExit(f"Transcript map has no transcript rows: {path}")

    if not saw_header:
        print(f"WARNING: No explicit transcript map header found in {path}; parsed as data rows.")
    return transcripts, proteins_present, genes_by_tx


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--transcript-map", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    catalog_ids, catalog_with_cds = catalog_transcripts(args.catalog)
    map_ids, with_protein, genes_by_tx = read_transcript_map(args.transcript_map)

    missing_rows = sorted(catalog_ids - map_ids)
    missing_protein_with_cds = sorted((catalog_ids - with_protein) & catalog_with_cds)
    missing_protein_no_cds = sorted((catalog_ids - with_protein) - catalog_with_cds)
    extra_rows = sorted(map_ids - catalog_ids)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", encoding="utf-8", newline="") as handle:
        handle.write("status\ttranscript_id\tgene\tdetail\n")
        for transcript in missing_rows:
            handle.write(f"missing_map_row\t{transcript}\t{genes_by_tx.get(transcript, '')}\t\n")
        for transcript in missing_protein_with_cds:
            handle.write(
                f"missing_protein_id_with_cds\t{transcript}\t{genes_by_tx.get(transcript, '')}\t"
                "catalog_has_cds\n"
            )
        for transcript in missing_protein_no_cds:
            handle.write(
                f"missing_protein_id_non_coding\t{transcript}\t{genes_by_tx.get(transcript, '')}\t"
                "catalog_has_no_cds\n"
            )
        for transcript in extra_rows:
            handle.write(f"extra_map_row\t{transcript}\t{genes_by_tx.get(transcript, '')}\t\n")

    failures = len(missing_rows) + len(missing_protein_with_cds)
    print(
        "{}: catalog_transcripts={} map_rows={} with_protein={} missing_rows={} missing_protein_with_cds={} missing_protein_no_cds={} extra_rows={}".format(
            args.report,
            len(catalog_ids),
            len(map_ids),
            len(with_protein),
            len(missing_rows),
            len(missing_protein_with_cds),
            len(missing_protein_no_cds),
            len(extra_rows),
        )
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())