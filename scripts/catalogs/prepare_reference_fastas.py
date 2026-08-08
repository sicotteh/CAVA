#!/usr/bin/env python3
"""Materialize gzip FASTA files and create .fai indexes without samtools.

The output is suitable for random-access validation by validate_gtf_sequences.py.
Coordinates in .fai are byte offsets into the uncompressed FASTA.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
from collections import defaultdict
from pathlib import Path


def materialize(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.stat().st_size > 0:
        print(f"Using existing {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    opener = gzip.open if source.suffix == ".gz" else open
    try:
        with opener(source, "rb") as src, temporary.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    temporary.replace(destination)
    print(destination)


def build_fai(fasta: Path) -> Path:
    index = Path(str(fasta) + ".fai")
    temporary = Path(str(index) + ".part")
    records: list[tuple[str, int, int, int, int]] = []
    with fasta.open("rb") as handle:
        name: str | None = None
        length = 0
        sequence_offset = 0
        line_bases = 0
        line_width = 0
        saw_short_line = False
        while True:
            line_offset = handle.tell()
            line = handle.readline()
            if not line:
                if name is not None:
                    records.append((name, length, sequence_offset, line_bases, line_width))
                break
            if line.startswith(b">"):
                if name is not None:
                    records.append((name, length, sequence_offset, line_bases, line_width))
                header = line[1:].strip().split(None, 1)[0]
                if not header:
                    raise SystemExit(f"Empty FASTA identifier at byte {line_offset}: {fasta}")
                name = header.decode("ascii", errors="strict")
                length = 0
                sequence_offset = handle.tell()
                line_bases = 0
                line_width = 0
                saw_short_line = False
                continue
            if name is None:
                raise SystemExit(f"Sequence data precedes first header: {fasta}")
            bases = len(line.rstrip(b"\r\n"))
            if bases == 0:
                continue
            if line_bases == 0:
                line_bases = bases
                line_width = len(line)
            elif bases != line_bases:
                saw_short_line = True
            elif saw_short_line:
                raise SystemExit(
                    f"Non-terminal short FASTA line for {name}; cannot create valid .fai: {fasta}"
                )
            length += bases
    if not records:
        raise SystemExit(f"No FASTA records found: {fasta}")
    with temporary.open("w", encoding="ascii", newline="\n") as out:
        for record in records:
            out.write("\t".join(str(item) for item in record) + "\n")
    temporary.replace(index)
    print(index)
    return index




def read_fasta_names(index: Path) -> set[str]:
    return {line.split("\t", 1)[0] for line in index.read_text(encoding="ascii").splitlines() if line}


def assembly_report_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9:
                rows.append(fields)
    return rows


def build_aliases(fasta: Path, reports: list[Path]) -> Path:
    index = Path(str(fasta) + ".fai")
    names = read_fasta_names(index)
    aliases: dict[str, set[str]] = defaultdict(set)
    for report in reports:
        for fields in assembly_report_rows(report):
            # Only use identifier-bearing assembly-report columns. Mapping
            # roles, units, lengths, or relationship fields as aliases could
            # resolve an unrelated GTF contig to the wrong FASTA record.
            identifier_indexes = (0, 4, 6, 9)
            values = {
                fields[index].strip()
                for index in identifier_indexes
                if index < len(fields)
                and fields[index].strip() not in {"na", "N/A", "-"}
            }
            matching_records = sorted(values & names)
            if len(matching_records) != 1:
                continue
            record = matching_records[0]
            for alias in values:
                aliases[alias].add(record)
                if alias.startswith("chr"):
                    aliases[alias[3:]].add(record)
                elif alias in {"X", "Y", "M", "MT"} or alias.isdigit():
                    aliases["chr" + ("M" if alias == "MT" else alias)].add(record)
    path = Path(str(fasta) + ".aliases.tsv")
    temporary = Path(str(path) + ".part")
    with temporary.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# alias\tfasta_record\n")
        for alias in sorted(aliases):
            records = aliases[alias]
            if len(records) == 1 and alias not in names:
                out.write("{}\t{}\n".format(alias, next(iter(records))))
    temporary.replace(path)
    print(path)
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, default=Path("catalog-build"))
    parser.add_argument("--source", type=Path, action="append", default=[])
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    source_dir = workspace / "fasta"
    output_dir = workspace / "prepared-fasta"
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.source:
        sources = [path.resolve() for path in args.source]
    else:
        sources = sorted(source_dir.glob("*.fa")) + sorted(source_dir.glob("*.fna"))
        sources += sorted(source_dir.glob("*.fa.gz")) + sorted(source_dir.glob("*.fna.gz"))
    if not sources:
        raise SystemExit(f"No FASTA files found under {source_dir}")
    reports = sorted(source_dir.glob("*_assembly_report.txt"))
    for source in sources:
        name = source.name[:-3] if source.suffix == ".gz" else source.name
        destination = output_dir / name
        materialize(source, destination)
        build_fai(destination)
        build_aliases(destination, reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
