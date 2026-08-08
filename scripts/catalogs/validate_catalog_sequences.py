#!/usr/bin/env python3
"""Validate a final CAVA catalog against genomic and transcript sequences.

CAVA catalog exon intervals are already 0-based, half-open. Minus-strand exon
sequences are reverse-complemented in transcript order before comparison with
the authoritative RNA FASTA.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from validate_gtf_sequences import FastaIndex, load_rna_sequences, reverse_complement


def read_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13 or (len(fields) - 11) % 2:
                raise SystemExit("Malformed catalog row {}:{}".format(path, number))
            yield number, fields


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--rna-fasta", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--require-gene", action="append", default=[])
    args = parser.parse_args()

    records = list(read_rows(args.catalog))
    identifiers = {fields[0] for _, fields in records}
    rna_exact, rna_by_base = load_rna_sequences(args.rna_fasta, identifiers)
    fasta = FastaIndex(args.genome_fasta)
    report_rows = []
    failures = 0
    found_genes: set[str] = set()
    try:
        for number, fields in records:
            transcript = fields[0]
            gene = fields[1]
            contig = fields[4]
            strand = fields[5]
            exons = [
                (int(fields[index]), int(fields[index + 1]))
                for index in range(11, len(fields), 2)
            ]
            status = "ok"
            detail = ""
            sequence_parts = []
            try:
                for start0, end0 in exons:
                    if start0 < 0 or end0 <= start0:
                        raise ValueError("invalid 0-based half-open exon {}-{}".format(start0, end0))
                    piece = fasta.fetch(contig, start0, end0)
                    sequence_parts.append(reverse_complement(piece) if strand == "-1" else piece)
            except (KeyError, ValueError) as exc:
                status = "genome_error"
                detail = str(exc)
            assembled = "".join(sequence_parts)
            rna_id = ""
            if status == "ok":
                if transcript in rna_exact:
                    rna_id = transcript
                else:
                    # load_rna_sequences indexes unique version/location-free IDs.
                    from validate_gtf_sequences import versionless
                    candidates = rna_by_base.get(versionless(transcript), [])
                    if len(candidates) == 1:
                        rna_id = candidates[0]
                    elif not candidates:
                        status = "missing_rna"
                    else:
                        status = "ambiguous_rna"
                        detail = ",".join(candidates)
                if rna_id and assembled != rna_exact[rna_id]:
                    status = "sequence_mismatch"
                    first = next(
                        (
                            index
                            for index, pair in enumerate(zip(assembled, rna_exact[rna_id]), 1)
                            if pair[0] != pair[1]
                        ),
                        0,
                    )
                    detail = "genomic_len={} rna_len={} first_mismatch_1based={}".format(
                        len(assembled), len(rna_exact[rna_id]), first
                    )
            if status == "ok":
                found_genes.add(gene.upper())
            else:
                failures += 1
            report_rows.append(
                {
                    "line": str(number),
                    "transcript_id": transcript,
                    "rna_id": rna_id,
                    "gene": gene,
                    "contig": contig,
                    "strand": strand,
                    "exon_count": str(len(exons)),
                    "status": status,
                    "detail": detail,
                }
            )
    finally:
        fasta.close()

    missing_genes = sorted({gene.upper() for gene in args.require_gene} - found_genes)
    for gene in missing_genes:
        failures += 1
        report_rows.append(
            {
                "line": "",
                "transcript_id": "",
                "rna_id": "",
                "gene": gene,
                "contig": "",
                "strand": "",
                "exon_count": "",
                "status": "missing_required_gene",
                "detail": "",
            }
        )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "line",
        "transcript_id",
        "rna_id",
        "gene",
        "contig",
        "strand",
        "exon_count",
        "status",
        "detail",
    ]
    with args.report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)
    print("{}: transcripts={} failures={}".format(args.report, len(records), failures))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
