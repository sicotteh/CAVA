#!/usr/bin/env python3
"""Create reviewable MANE 1.5 GRCh37 candidate transcript lists.

RefSeq candidates are chosen from NCBI's GRCh37 annotation, preferring an exact
accession.version and otherwise the highest earlier version of the same NM_/NR_
accession. Ensembl candidates are matched by gene and ordered CDS exon lengths,
then require exact coding-sequence identity across GRCh38 and GRCh37. UTR
lengths are intentionally ignored.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from validate_gtf_sequences import FastaIndex, assembled_sequence, open_text, parse_attributes, Transcript


def split_version(identifier: str) -> tuple[str, int | None]:
    match = re.match(r"^(.*)\.(\d+)$", identifier)
    return (match.group(1), int(match.group(2))) if match else (identifier, None)


def read_mane_summary(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        rows = [{key: (value or "") for key, value in row.items()} for row in reader]
    if not rows:
        raise SystemExit(f"No MANE rows found: {path}")
    return rows


def find_column(columns: list[str], patterns: tuple[str, ...]) -> str:
    lowered = {column.lower(): column for column in columns}
    for pattern in patterns:
        for lower, original in lowered.items():
            if pattern in lower:
                return original
    raise SystemExit(f"Could not identify MANE summary column matching {patterns}; columns={columns}")


def parse_gtf_transcripts(path: Path) -> dict[str, Transcript]:
    transcripts: dict[str, Transcript] = {}
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            attrs = parse_attributes(fields[8])
            identifier = attrs.get("transcript_id")
            if not identifier:
                continue
            tx = transcripts.setdefault(
                identifier,
                Transcript(identifier, fields[0], fields[6], attrs.get("gene_name", attrs.get("gene", ""))),
            )
            interval = (int(fields[3]) - 1, int(fields[4]))
            if fields[2] == "exon":
                tx.exons.add(interval)
            elif fields[2] == "CDS":
                tx.cds.add(interval)
    return transcripts


def cds_signature(tx: Transcript) -> tuple[int, ...]:
    return tuple(end - start for start, end in tx.intervals(cds_only=True))


def choose_refseq(target: str, available: dict[str, list[tuple[int | None, str]]]) -> tuple[str, str]:
    base, version = split_version(target)
    candidates = available.get(base, [])
    if not candidates:
        return "", "no_accession_on_grch37"
    exact = [identifier for candidate_version, identifier in candidates if candidate_version == version]
    if exact:
        return exact[0], "exact_accession_version"
    earlier = [(candidate_version, identifier) for candidate_version, identifier in candidates if candidate_version is not None and (version is None or candidate_version < version)]
    if earlier:
        return max(earlier)[1], "earlier_accession_version"
    return "", "only_newer_or_unversioned_candidates"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mane-summary", type=Path, required=True)
    parser.add_argument("--mane-refseq-grch38-gtf", type=Path, required=True)
    parser.add_argument("--mane-ensembl-grch38-gtf", type=Path, required=True)
    parser.add_argument("--ensembl75-grch37-gtf", type=Path, required=True)
    parser.add_argument("--grch38-refseq-fasta", type=Path, required=True)
    parser.add_argument("--grch38-fasta", type=Path, required=True, help="GRCh38 FASTA matching the MANE Ensembl GTF")
    parser.add_argument("--grch37-fasta", type=Path, required=True)
    parser.add_argument("--grch37-refseq-gtf", type=Path, required=True)
    parser.add_argument("--refseq-output", type=Path, required=True)
    parser.add_argument("--ensembl-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    summary = read_mane_summary(args.mane_summary)
    columns = list(summary[0])
    refseq_col = find_column(columns, ("refseq_nuc", "refseq nucleotide", "refseq"))
    ensembl_col = find_column(columns, ("ensembl_nuc", "ensembl nucleotide", "ensembl"))
    symbol_col = find_column(columns, ("symbol", "gene name", "gene"))

    gr38_refseq = parse_gtf_transcripts(args.mane_refseq_grch38_gtf)
    gr37_refseq = parse_gtf_transcripts(args.grch37_refseq_gtf)
    refseq_by_base: dict[str, list[tuple[int | None, str]]] = defaultdict(list)
    for identifier in gr37_refseq:
        if identifier.startswith(("NM_", "NR_", "XM_", "XR_")):
            base, version = split_version(identifier)
            refseq_by_base[base].append((version, identifier))

    gr38_ensembl = parse_gtf_transcripts(args.mane_ensembl_grch38_gtf)
    gr37_ensembl = parse_gtf_transcripts(args.ensembl75_grch37_gtf)
    gr37_by_gene_signature: dict[tuple[str, tuple[int, ...]], list[Transcript]] = defaultdict(list)
    for tx in gr37_ensembl.values():
        if tx.cds and tx.gene_name:
            gr37_by_gene_signature[(tx.gene_name.upper(), cds_signature(tx))].append(tx)

    fasta38_refseq = FastaIndex(args.grch38_refseq_fasta)
    fasta38 = FastaIndex(args.grch38_fasta)
    fasta37 = FastaIndex(args.grch37_fasta)
    report_rows: list[dict[str, str]] = []
    accepted_refseq: set[str] = set()
    accepted_ensembl: set[str] = set()
    try:
        for row in summary:
            target_refseq = row.get(refseq_col, "").strip()
            target_ensembl = row.get(ensembl_col, "").strip()
            symbol = row.get(symbol_col, "").strip()
            refseq_choice, refseq_method = choose_refseq(target_refseq, refseq_by_base) if target_refseq else ("", "missing_mane_refseq")
            refseq_signature_match = ""
            refseq_sequence_match = ""
            if refseq_choice:
                tx38_refseq = gr38_refseq.get(target_refseq)
                tx37_refseq = gr37_refseq.get(refseq_choice)
                if tx38_refseq is None:
                    refseq_method += ":grch38_definition_missing"
                    refseq_choice = ""
                elif tx37_refseq is None:
                    refseq_method += ":grch37_definition_missing"
                    refseq_choice = ""
                elif not tx38_refseq.cds or not tx37_refseq.cds:
                    refseq_method += ":noncoding_or_missing_cds_not_auto_mapped"
                    refseq_choice = ""
                else:
                    refseq_signature_match = str(cds_signature(tx38_refseq) == cds_signature(tx37_refseq)).lower()
                    if refseq_signature_match != "true":
                        refseq_method += ":coding_exon_lengths_differ"
                        refseq_choice = ""
                    else:
                        try:
                            sequence38_refseq = assembled_sequence(tx38_refseq, fasta38_refseq, cds_only=True)
                            sequence37_refseq = assembled_sequence(tx37_refseq, fasta37, cds_only=True)
                        except (KeyError, ValueError):
                            refseq_method += ":cds_sequence_unavailable"
                            refseq_choice = ""
                        else:
                            refseq_sequence_match = str(sequence38_refseq == sequence37_refseq).lower()
                            if refseq_sequence_match == "true":
                                accepted_refseq.add(refseq_choice)
                                refseq_method += ":cds_lengths_and_sequence_match"
                            else:
                                refseq_method += ":cds_sequence_differs"
                                refseq_choice = ""

            ensembl_choice = ""
            ensembl_method = ""
            candidates: list[Transcript] = []
            tx38 = gr38_ensembl.get(target_ensembl)
            if tx38 is None:
                ensembl_method = "grch38_transcript_missing_from_gtf"
            elif not tx38.cds:
                ensembl_method = "noncoding_or_missing_cds_not_auto_mapped"
            else:
                gene = (tx38.gene_name or symbol).upper()
                candidates = gr37_by_gene_signature.get((gene, cds_signature(tx38)), [])
                try:
                    sequence38 = assembled_sequence(tx38, fasta38, cds_only=True)
                except (KeyError, ValueError):
                    sequence38 = ""
                    ensembl_method = "grch38_cds_sequence_unavailable"
                exact: list[str] = []
                if sequence38:
                    for candidate in candidates:
                        try:
                            if assembled_sequence(candidate, fasta37, cds_only=True) == sequence38:
                                exact.append(candidate.identifier)
                        except (KeyError, ValueError):
                            continue
                    if len(exact) == 1:
                        ensembl_choice = exact[0]
                        accepted_ensembl.add(ensembl_choice)
                        ensembl_method = "unique_gene_cds_lengths_and_exact_cds_sequence"
                    elif len(exact) > 1:
                        ensembl_method = "ambiguous_exact_cds_matches:" + ",".join(sorted(exact))
                    elif not ensembl_method:
                        ensembl_method = "no_exact_cds_sequence_match"

            report_rows.append(
                {
                    "gene": symbol,
                    "mane_refseq_grch38": target_refseq,
                    "refseq_grch37_candidate": refseq_choice,
                    "refseq_method": refseq_method,
                    "refseq_cds_length_signature_match": refseq_signature_match,
                    "refseq_exact_cds_sequence_match": refseq_sequence_match,
                    "mane_ensembl_grch38": target_ensembl,
                    "ensembl_grch37_candidate": ensembl_choice,
                    "ensembl_method": ensembl_method,
                    "ensembl_length_signature_candidates": ",".join(sorted(tx.identifier for tx in candidates)),
                }
            )
    finally:
        fasta38_refseq.close()
        fasta38.close()
        fasta37.close()

    for path, identifiers in ((args.refseq_output, accepted_refseq), (args.ensembl_output, accepted_ensembl)):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(identifier + "\n" for identifier in sorted(identifiers)), encoding="utf-8")
        print(f"{path}: {len(identifiers)} candidates")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(report_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)
    print(args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
