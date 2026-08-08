#!/usr/bin/env python3
"""Validate and select transcript loci using genomic and RNA sequences.

GTF coordinates are 1-based inclusive. FASTA slicing is 0-based half-open, so
GTF intervals are converted to ``[start - 1, end)``. Exons are joined in
transcript order and minus-strand pieces are reverse-complemented.

When a transcript accession is annotated at multiple genomic loci, every locus
is checked independently. One exact RNA-sequence match is selected
deterministically, preferring a primary contig. ``--filtered-gtf`` writes only
the selected locus, preventing CAVA's builder from conflating duplicate
transcript IDs across reference and alternate loci.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ATTR_RE = re.compile(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"')
RC = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn"
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC)[::-1]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def canonical_identifier(identifier: str) -> str:
    """Remove a MANE location suffix and then an accession version."""
    identifier = re.sub(r"_[0-9]+$", "", identifier)
    return identifier.rsplit(".", 1)[0] if re.search(r"\.\d+$", identifier) else identifier


# Kept as a public helper for the other catalog scripts.
versionless = canonical_identifier


class FastaIndex:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        index_path = Path(str(fasta) + ".fai")
        if not index_path.is_file():
            raise SystemExit(
                "Missing FASTA index {}; run prepare_reference_fastas.py".format(
                    index_path
                )
            )
        self.records: dict[str, tuple[int, int, int, int]] = {}
        for line in index_path.read_text(encoding="ascii").splitlines():
            name, length, offset, line_bases, line_width = line.split("\t")[:5]
            self.records[name] = tuple(
                map(int, (length, offset, line_bases, line_width))
            )
        self.aliases: dict[str, set[str]] = defaultdict(set)
        for name in self.records:
            self._add_alias(name, name)
            if name.startswith("chr"):
                self._add_alias(name[3:], name)
            else:
                self._add_alias("chr" + name, name)
        self._read_header_aliases()
        self._read_sidecar_aliases()
        self.handle = fasta.open("rb")

    def _add_alias(self, alias: str, name: str) -> None:
        if alias:
            self.aliases[alias].add(name)

    def _read_header_aliases(self) -> None:
        chromosome_re = re.compile(rb"\bchromosome\s+([0-9XYM]+)\b", re.I)
        with self.fasta.open("rb") as handle:
            for line in handle:
                if not line.startswith(b">"):
                    continue
                parts = line[1:].rstrip().split(None, 1)
                name = parts[0].decode("ascii", errors="ignore")
                description = parts[1] if len(parts) > 1 else b""
                lowered = description.lower()
                # Only add chr-style aliases from primary assembled chromosomes.
                # Unlocalized, unplaced, alternate, fix, patch, and scaffold
                # records also mention chromosome numbers and would otherwise make
                # aliases like chr1 ambiguous.
                if b"primary assembly" not in lowered:
                    continue
                if any(
                    token in lowered
                    for token in (
                        b"unlocalized",
                        b"unplaced",
                        b"alternate",
                        b"scaffold",
                        b"patch",
                    )
                ):
                    continue
                match = chromosome_re.search(description)
                if match:
                    chromosome = match.group(1).decode("ascii").upper()
                    if chromosome == "MT":
                        chromosome = "M"
                    self._add_alias(chromosome, name)
                    self._add_alias("chr" + chromosome, name)

    def _read_sidecar_aliases(self) -> None:
        alias_path = Path(str(self.fasta) + ".aliases.tsv")
        if not alias_path.is_file():
            return
        for number, line in enumerate(
            alias_path.read_text(encoding="utf-8").splitlines(), 1
        ):
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise SystemExit(
                    "Malformed FASTA alias at {}:{}".format(alias_path, number)
                )
            alias, record = fields
            if record not in self.records:
                raise SystemExit(
                    "FASTA alias points to an unknown record at {}:{}: {}".format(
                        alias_path, number, record
                    )
                )
            self._add_alias(alias, record)

    def resolve(self, requested: str) -> str:
        if requested in self.records:
            return requested
        candidates = self.aliases.get(requested, set())
        if len(candidates) == 1:
            return next(iter(candidates))
        raise KeyError(requested)

    def fetch(self, contig: str, start0: int, end0: int) -> str:
        name = self.resolve(contig)
        length, offset, line_bases, line_width = self.records[name]
        if not (0 <= start0 <= end0 <= length):
            raise ValueError(
                "Invalid interval {}:{}-{}; length={}".format(
                    contig, start0, end0, length
                )
            )
        pieces: list[bytes] = []
        position = start0
        while position < end0:
            line_number, in_line = divmod(position, line_bases)
            take = min(end0 - position, line_bases - in_line)
            self.handle.seek(offset + line_number * line_width + in_line)
            pieces.append(self.handle.read(take))
            position += take
        return b"".join(pieces).decode("ascii").upper()

    def close(self) -> None:
        self.handle.close()


@dataclass
class Transcript:
    identifier: str
    contig: str
    strand: str
    gene_name: str = ""
    exons: set[tuple[int, int]] = field(default_factory=set)
    cds: set[tuple[int, int]] = field(default_factory=set)

    def intervals(self, cds_only: bool = False) -> list[tuple[int, int]]:
        values = self.cds if cds_only else self.exons
        return sorted(values, reverse=self.strand == "-")


def parse_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def read_transcript_filter(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }


def parse_gtf(
    path: Path, selected: set[str] | None
) -> dict[str, dict[tuple[str, str], Transcript]]:
    result: dict[str, dict[tuple[str, str], Transcript]] = defaultdict(dict)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            attrs = parse_attributes(fields[8])
            identifier = attrs.get("transcript_id")
            if not identifier or (selected is not None and identifier not in selected):
                continue
            start1, end1 = int(fields[3]), int(fields[4])
            if start1 < 1 or end1 < start1:
                raise SystemExit(
                    "Invalid 1-based GTF interval at {}:{}".format(path, line_number)
                )
            key = (fields[0], fields[6])
            transcript = result[identifier].setdefault(
                key,
                Transcript(
                    identifier,
                    fields[0],
                    fields[6],
                    attrs.get("gene_name", attrs.get("gene", attrs.get("gene_id", ""))),
                ),
            )
            interval = (start1 - 1, end1)
            if fields[2] == "exon":
                transcript.exons.add(interval)
            elif fields[2] == "CDS":
                transcript.cds.add(interval)
    return dict(result)


def load_rna_sequences(
    path: Path, selected: set[str]
) -> tuple[dict[str, str], dict[str, list[str]]]:
    exact: dict[str, str] = {}
    by_base: dict[str, list[str]] = defaultdict(list)
    current: str | None = None
    chunks: list[str] = []

    def save() -> None:
        if current is None:
            return
        exact[current] = "".join(chunks).upper().replace("U", "T")
        by_base[canonical_identifier(current)].append(current)

    selected_bases = {canonical_identifier(item) for item in selected}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                tokens = re.split(r"[|\s]", line[1:].strip())
                current = next(
                    (
                        token
                        for token in tokens
                        if token in selected
                        or canonical_identifier(token) in selected_bases
                    ),
                    None,
                )
                chunks = []
            elif current is not None:
                chunks.append(line.strip())
        save()
    return exact, by_base


def assembled_sequence(
    transcript: Transcript, fasta: FastaIndex, cds_only: bool = False
) -> str:
    chunks: list[str] = []
    for start0, end0 in transcript.intervals(cds_only=cds_only):
        sequence = fasta.fetch(transcript.contig, start0, end0)
        chunks.append(
            reverse_complement(sequence) if transcript.strand == "-" else sequence
        )
    return "".join(chunks)


def primary_contig_rank(contig: str) -> int:
    value = contig[3:] if contig.startswith("chr") else contig
    if "_" not in value and value in {
        *(str(number) for number in range(1, 23)),
        "X",
        "Y",
        "M",
        "MT",
    }:
        return 0
    if re.fullmatch(r"NC_0*(?:[1-9]|1[0-9]|2[0-4])\.\d+", contig):
        return 0
    if contig.startswith("NC_012920."):
        return 0
    return 1


def resolve_rna_id(
    identifier: str, exact: dict[str, str], by_base: dict[str, list[str]]
) -> tuple[str, str]:
    if identifier in exact:
        return identifier, ""
    candidates = by_base.get(canonical_identifier(identifier), [])
    if len(candidates) == 1:
        return candidates[0], ""
    if not candidates:
        return "", "No exact or unique version/location-free RNA sequence"
    return "", "Ambiguous RNA sequences: " + ",".join(sorted(candidates))


def write_filtered_gtf(
    source: Path,
    destination: Path,
    selected_loci: dict[str, tuple[str, str]],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if destination.suffix == ".gz" else open
    with open_text(source) as handle, opener(
        destination, "wt", encoding="utf-8", newline="\n"
    ) as out:
        for line in handle:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attributes(fields[8])
            identifier = attrs.get("transcript_id")
            if identifier and selected_loci.get(identifier) == (fields[0], fields[6]):
                out.write(line)
    print(destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--rna-fasta", type=Path)
    parser.add_argument("--transcripts", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--allow-mismatches", action="store_true")
    parser.add_argument("--require-gene", action="append", default=[])
    parser.add_argument("--max-transcripts", type=int, default=0)
    parser.add_argument("--accepted-transcripts", type=Path)
    parser.add_argument("--filtered-gtf", type=Path)
    args = parser.parse_args()

    selected = read_transcript_filter(args.transcripts)
    transcript_loci = parse_gtf(args.gtf, selected)
    missing_definitions = sorted((selected or set()) - set(transcript_loci))
    ordered_ids = sorted(transcript_loci)
    if args.max_transcripts:
        ordered_ids = ordered_ids[: args.max_transcripts]

    rna_exact: dict[str, str] = {}
    rna_by_base: dict[str, list[str]] = {}
    if args.rna_fasta:
        rna_exact, rna_by_base = load_rna_sequences(args.rna_fasta, set(ordered_ids))

    fasta = FastaIndex(args.genome_fasta)
    rows: list[dict[str, str]] = []
    failed_transcripts = 0
    accepted_ids: list[str] = []
    selected_loci: dict[str, tuple[str, str]] = {}
    accepted_genes: set[str] = set()
    try:
        for identifier in ordered_ids:
            rna_id, rna_error = (
                resolve_rna_id(identifier, rna_exact, rna_by_base)
                if args.rna_fasta
                else ("", "")
            )
            candidates = []
            for locus_key, transcript in transcript_loci[identifier].items():
                status = "ok"
                detail = ""
                assembled = ""
                try:
                    if not transcript.exons:
                        status, detail = "missing_exons", "No exon records"
                    else:
                        assembled = assembled_sequence(transcript, fasta)
                except KeyError:
                    status, detail = "missing_contig", transcript.contig
                except ValueError as exc:
                    status, detail = "invalid_interval", str(exc)

                if status == "ok" and args.rna_fasta:
                    if not rna_id:
                        status, detail = "missing_or_ambiguous_rna", rna_error
                    elif assembled != rna_exact[rna_id]:
                        status = "sequence_mismatch"
                        first = next(
                            (
                                index
                                for index, pair in enumerate(
                                    zip(assembled, rna_exact[rna_id]), 1
                                )
                                if pair[0] != pair[1]
                            ),
                            0,
                        )
                        detail = (
                            "genomic_len={} rna_len={} first_mismatch_1based={}".format(
                                len(assembled), len(rna_exact[rna_id]), first
                            )
                        )
                candidates.append((locus_key, transcript, status, detail, rna_id))

            exact_candidates = [item for item in candidates if item[2] == "ok"]
            chosen = None
            if exact_candidates:
                chosen = min(
                    exact_candidates,
                    key=lambda item: (
                        primary_contig_rank(item[1].contig),
                        item[1].contig,
                        item[1].strand,
                    ),
                )
                selected_loci[identifier] = chosen[0]
                accepted_ids.append(identifier)
                if chosen[1].gene_name:
                    accepted_genes.add(chosen[1].gene_name.upper())
            else:
                failed_transcripts += 1

            for locus_key, transcript, status, detail, candidate_rna_id in candidates:
                selected_flag = chosen is not None and locus_key == chosen[0]
                reported_status = status
                if status == "ok" and not selected_flag:
                    reported_status = "exact_duplicate_not_selected"
                    detail = "A deterministic preferred exact locus was selected"
                rows.append(
                    {
                        "transcript_id": identifier,
                        "rna_id": candidate_rna_id,
                        "gene_name": transcript.gene_name,
                        "contig": transcript.contig,
                        "strand": transcript.strand,
                        "exon_count": str(len(transcript.exons)),
                        "cds_exon_count": str(len(transcript.cds)),
                        "selected": "yes" if selected_flag else "no",
                        "status": reported_status,
                        "detail": detail,
                    }
                )
    finally:
        fasta.close()

    for identifier in missing_definitions:
        failed_transcripts += 1
        rows.append(
            {
                "transcript_id": identifier,
                "rna_id": "",
                "gene_name": "",
                "contig": "",
                "strand": "",
                "exon_count": "0",
                "cds_exon_count": "0",
                "selected": "no",
                "status": "missing_gtf",
                "detail": "Transcript list entry is absent from GTF",
            }
        )

    required_genes = {gene.upper() for gene in args.require_gene}
    missing_required_genes = sorted(required_genes - accepted_genes)
    for gene in missing_required_genes:
        rows.append(
            {
                "transcript_id": "",
                "rna_id": "",
                "gene_name": gene,
                "contig": "",
                "strand": "",
                "exon_count": "",
                "cds_exon_count": "",
                "selected": "no",
                "status": "missing_required_gene",
                "detail": "Required gene has no exact reference-consistent locus",
            }
        )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "transcript_id",
        "rna_id",
        "gene_name",
        "contig",
        "strand",
        "exon_count",
        "cds_exon_count",
        "selected",
        "status",
        "detail",
    ]
    with args.report.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if args.accepted_transcripts:
        args.accepted_transcripts.parent.mkdir(parents=True, exist_ok=True)
        args.accepted_transcripts.write_text(
            "".join(identifier + "\n" for identifier in accepted_ids),
            encoding="utf-8",
        )
        print(
            "{}: accepted={}".format(args.accepted_transcripts, len(accepted_ids))
        )
    if args.filtered_gtf:
        write_filtered_gtf(args.gtf, args.filtered_gtf, selected_loci)

    print(
        "{}: transcripts={} accepted={} failed={}".format(
            args.report, len(ordered_ids), len(accepted_ids), failed_transcripts
        )
    )
    if missing_required_genes:
        return 1
    if failed_transcripts and not args.allow_mismatches:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
