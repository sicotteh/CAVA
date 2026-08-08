#!/usr/bin/env python3
"""Keep only catalog transcripts with non-empty protein IDs in companion map."""
from __future__ import annotations

import argparse
import gzip
import tempfile
from pathlib import Path

import pysam


def read_map(path: Path):
    comments: list[str] = []
    header: str | None = None
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                continue
            fields = line.split("\t")
            if header is None:
                header = line
                if len(fields) < 4 or fields[:4] != ["GENEID", "SYMBOL", "Transcript", "Protein"]:
                    raise SystemExit(f"Unexpected transcript map header in {path}: {fields}")
                continue
            if len(fields) < 4:
                raise SystemExit(f"Malformed transcript map row in {path}: {fields}")
            rows.append(fields)
    if header is None:
        raise SystemExit(f"Transcript map is empty: {path}")
    return comments, header, rows


def read_catalog_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header_lines: list[str] = []
    body: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                body.append(line)
    return header_lines, body


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--transcript-map", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    comments, header, map_rows = read_map(args.transcript_map)
    with_protein = {row[2].strip() for row in map_rows if row[3].strip()}

    catalog_headers, catalog_rows = read_catalog_rows(args.catalog)
    retained_catalog: list[str] = []
    dropped_catalog = 0
    for line in catalog_rows:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 13:
            raise SystemExit(f"Malformed catalog row in {args.catalog}: {line[:200]}")
        if fields[0] in with_protein:
            retained_catalog.append(line)
        else:
            dropped_catalog += 1

    retained_map = [row for row in map_rows if row[2].strip() in with_protein]
    dropped_map = len(map_rows) - len(retained_map)

    with tempfile.TemporaryDirectory(prefix="cava-enforce-protein-") as tmp_dir:
        tmp = Path(tmp_dir)

        # Rewrite transcript map in place.
        map_tmp = tmp / "map.txt"
        with map_tmp.open("w", encoding="utf-8", newline="\n") as out:
            for line in comments:
                out.write(line + "\n")
            out.write(header + "\n")
            for row in retained_map:
                out.write("\t".join(row) + "\n")
        map_tmp.replace(args.transcript_map)

        # Rewrite compressed catalog in place and rebuild index.
        plain = tmp / "catalog.tsv"
        with plain.open("w", encoding="utf-8", newline="\n") as out:
            for line in catalog_headers:
                out.write(line)
            for line in retained_catalog:
                out.write(line)

        catalog_path = args.catalog
        index_path = Path(str(catalog_path) + ".tbi")
        if catalog_path.suffix == ".gz":
            pysam.tabix_compress(str(plain), str(catalog_path), force=True)
        else:
            plain.replace(catalog_path)
        pysam.tabix_index(
            str(catalog_path),
            seq_col=4,
            start_col=6,
            end_col=7,
            zerobased=True,
            meta_char="#",
            force=True,
        )
        if not index_path.is_file() or index_path.stat().st_size == 0:
            raise SystemExit(f"Failed to rebuild index: {index_path}")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        "\n".join(
            [
                "metric\tvalue",
                f"catalog_rows_total\t{len(catalog_rows)}",
                f"catalog_rows_retained\t{len(retained_catalog)}",
                f"catalog_rows_dropped_missing_protein\t{dropped_catalog}",
                f"map_rows_total\t{len(map_rows)}",
                f"map_rows_retained\t{len(retained_map)}",
                f"map_rows_dropped_missing_protein\t{dropped_map}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        f"{args.report}: catalog_retained={len(retained_catalog)}/{len(catalog_rows)} "
        f"map_retained={len(retained_map)}/{len(map_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())