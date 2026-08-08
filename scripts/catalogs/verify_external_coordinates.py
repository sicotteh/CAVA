#!/usr/bin/env python3
"""Cross-check CAVA catalog transcript coordinates against Ensembl REST.

Use this validator only when external confirmation is warranted:
- during active build/debug phases,
- after coordinate parsing/transform code changes,
- or when runtime discrepancies suggest coordinate issues.

This catches incorrect transforms that can pass internal consistency checks.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import subprocess
from pathlib import Path


def versionless(transcript_id: str) -> str:
    return transcript_id.split(".", 1)[0]


def read_catalog(path: Path) -> dict[str, dict[str, str | int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    rows: dict[str, dict[str, str | int]] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise SystemExit(f"Malformed catalog row at {path}:{line_number}")
            tid = fields[0]
            base = versionless(tid)
            if base in rows:
                continue
            rows[base] = {
                "transcript": tid,
                "gene": fields[1],
                "chrom": fields[4],
                "strand": fields[5],
                "start0": int(fields[6]),
                "end1": int(fields[7]),
            }
    if not rows:
        raise SystemExit(f"No transcript rows found in {path}")
    return rows


def lookup_ensembl(ensembl_id: str, curl_bin: str, timeout_seconds: int) -> dict[str, object]:
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    result = subprocess.run(
        [curl_bin, "-sS", "-L", "--max-time", str(timeout_seconds), url],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"curl exit {result.returncode}")
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError("Invalid JSON from Ensembl REST") from exc
    if "error" in payload:
        raise RuntimeError(str(payload["error"]))
    return payload


def pick_transcripts(
    catalog_rows: dict[str, dict[str, str | int]],
    requested: list[str],
    sample_size: int,
    seed: int,
) -> list[str]:
    selected: list[str] = []
    seen = set()
    for item in requested:
        base = versionless(item)
        if base in catalog_rows and base not in seen:
            selected.append(base)
            seen.add(base)
    remaining = [key for key in sorted(catalog_rows) if key not in seen]
    if sample_size > 0 and remaining:
        random.seed(seed)
        k = min(sample_size, len(remaining))
        selected.extend(sorted(random.sample(remaining, k)))
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=25)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--transcript", action="append", default=[])
    parser.add_argument("--curl-bin", default="curl")
    parser.add_argument("--timeout-seconds", type=int, default=30)
    args = parser.parse_args()

    rows = read_catalog(args.catalog)
    selected = pick_transcripts(rows, args.transcript, args.sample_size, args.seed)
    if not selected:
        raise SystemExit("No transcripts were selected for external verification")

    checked = 0
    failures = 0
    output_rows: list[dict[str, str]] = []

    for base_id in selected:
        local = rows[base_id]
        status = "ok"
        detail = ""
        rest_chrom = ""
        rest_start = ""
        rest_end = ""
        rest_strand = ""
        try:
            payload = lookup_ensembl(base_id, args.curl_bin, args.timeout_seconds)
            rest_chrom = str(payload.get("seq_region_name", ""))
            rest_start = str(payload.get("start", ""))
            rest_end = str(payload.get("end", ""))
            rest_strand = str(payload.get("strand", ""))

            expected_start1 = int(local["start0"]) + 1
            expected_end1 = int(local["end1"])
            local_strand = int(local["strand"])
            api_start = int(rest_start)
            api_end = int(rest_end)
            api_strand = int(rest_strand)

            coord_ok = expected_start1 == api_start and expected_end1 == api_end
            strand_ok = local_strand == api_strand
            if coord_ok and strand_ok:
                status = "ok"
            elif not coord_ok and strand_ok:
                status = "coordinate_mismatch"
                detail = (
                    f"catalog_start1={expected_start1},rest_start1={api_start},"
                    f"catalog_end1={expected_end1},rest_end1={api_end}"
                )
            elif coord_ok and not strand_ok:
                status = "strand_mismatch"
                detail = f"catalog_strand={local_strand},rest_strand={api_strand}"
            else:
                status = "coordinate_and_strand_mismatch"
                detail = (
                    f"catalog_start1={expected_start1},rest_start1={api_start},"
                    f"catalog_end1={expected_end1},rest_end1={api_end},"
                    f"catalog_strand={local_strand},rest_strand={api_strand}"
                )
        except Exception as exc:
            status = "lookup_error"
            detail = str(exc)

        checked += 1
        if status != "ok":
            failures += 1

        output_rows.append(
            {
                "status": status,
                "transcript_id": str(local["transcript"]),
                "gene": str(local["gene"]),
                "catalog_chrom": str(local["chrom"]),
                "catalog_strand": str(local["strand"]),
                "catalog_start0": str(local["start0"]),
                "catalog_end1": str(local["end1"]),
                "ensembl_chrom": rest_chrom,
                "ensembl_strand": rest_strand,
                "ensembl_start1": rest_start,
                "ensembl_end1": rest_end,
                "detail": detail,
            }
        )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "status",
        "transcript_id",
        "gene",
        "catalog_chrom",
        "catalog_strand",
        "catalog_start0",
        "catalog_end1",
        "ensembl_chrom",
        "ensembl_strand",
        "ensembl_start1",
        "ensembl_end1",
        "detail",
    ]
    with args.report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"{args.report}: checked={checked} failures={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
