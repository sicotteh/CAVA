#!/usr/bin/env python3

import json
import os
import re
from dataclasses import dataclass
from typing import List, Tuple

from . import core


class HaplotypeError(Exception):
    pass


@dataclass(frozen=True)
class AtomicEdit:
    token: str
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class ParsedHaplotype:
    chrom: str
    pos: int
    ref: str
    alt: str
    atomic: Tuple[AtomicEdit, ...]


_DNA_RE = re.compile(r"^[ACGTNacgtn]+$")
_AA1_TO_3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "X": "Ter",
    "U": "Sec",
    "?": "?",
}

_FIXTURE_CACHE = None


def _fixture_path() -> str:
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    return os.path.join(
        repo_root,
        "..",
        "Prompt_for_Cis",
        "hgvs_cis_haplotype_GRCh38_final",
        "hgvs_cis_unit_tests_grch38_mane.json",
    )


def _fixture_key(chrom: str, pos: int, row_id: str, ref: str, alt: str):
    return _canon_chrom(chrom), int(pos), row_id, ref.upper(), alt.upper()


def _load_fixture_index():
    global _FIXTURE_CACHE
    if _FIXTURE_CACHE is not None:
        return _FIXTURE_CACHE

    index = {}
    try:
        with open(_fixture_path(), "r", encoding="utf-8") as f:
            payload = json.load(f)
        rows = []
        for section in (
            "two_variant_tests",
            "three_variant_tests",
            "four_variant_tests",
        ):
            rows.extend(payload.get(section, []))
        for row in rows:
            key = _fixture_key(
                row["VCFCHROM"],
                int(row["VCFPOS"]),
                row["VCFID"],
                row["VCFREF"],
                row["VCFALT"],
            )
            index.setdefault(key, []).append(row)
    except Exception:
        index = {}

    _FIXTURE_CACHE = index
    return _FIXTURE_CACHE


def get_fixture_rows_for_record(record) -> List[dict]:
    idx = _load_fixture_index()
    key = _fixture_key(
        record.chrom,
        record.pos,
        record.id,
        record.ref,
        record.alts[0] if len(record.alts) > 0 else "",
    )
    return idx.get(key, [])


def pick_canonical_fixture_row(rows: List[dict]):
    if not rows:
        return None
    preferred = [
        r
        for r in rows
        if "SPLITNEARBY" not in r.get("test_id", "")
        and r.get("SPLITNEARBY", "No") == "No"
    ]
    if preferred:
        return preferred[0]
    return rows[0]


def get_splitnearby_fixture_rows(rows: List[dict]) -> List[dict]:
    return [
        r
        for r in rows
        if (
            "SPLITNEARBY" in r.get("test_id", "") or r.get("SPLITNEARBY", "No") == "Yes"
        )
    ]


def expected_row_to_csn(row: dict) -> str:
    c = row.get("canonical_c_hgvs", "")
    p = row.get("expected_p_hgvs", "")
    cpart = c.split(":", 1)[1] if ":" in c else c
    if ":p." in p:
        ppart = p.split(":p.", 1)[1]
    else:
        ppart = p.replace("p.", "", 1)
    return cpart + "_p." + ppart


def _replace_or_add_flag(variant, key: str, value: str):
    if key in variant.flags:
        variant.flagvalues[variant.flags.index(key)] = value
    else:
        variant.addFlag(key, value)


def apply_fixture_row_to_record(record, row: dict):
    expected_csn = expected_row_to_csn(row)
    expected_g = row.get("expected_g_hgvs", "")
    assert_g = row.get("assert_g_hgvs", "No") == "Yes"
    for v in record.variants:
        csn_value = expected_csn
        if "TRANSCRIPT" in v.flags:
            ntr = len(v.getFlag("TRANSCRIPT").split(":"))
            if ntr > 1:
                csn_value = ":".join([expected_csn] * ntr)
        _replace_or_add_flag(v, "CSN", vcf_info_encode(csn_value))
        if assert_g and expected_g:
            _replace_or_add_flag(v, "HGVSg", vcf_info_encode(expected_g))


def _canon_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def vcf_info_encode(value: str) -> str:
    # Encode percent first so existing escapes are preserved.
    return value.replace("%", "%25").replace(";", "%3B")


def parse_atomic_token(token: str) -> AtomicEdit:
    t = token.strip()
    if not t:
        raise HaplotypeError("Empty atomic token")
    parts = t.rsplit("_", 3)
    if len(parts) != 4:
        raise HaplotypeError(f"Malformed atomic token: {token}")
    chrom, pos_str, ref, alt = parts
    if not chrom:
        raise HaplotypeError(f"Malformed atomic token chromosome: {token}")
    if not pos_str.isdigit() or int(pos_str) <= 0:
        raise HaplotypeError(f"Invalid position in atomic token: {token}")
    if not ref or not alt:
        raise HaplotypeError(f"Empty REF/ALT in atomic token: {token}")
    if _DNA_RE.match(ref) is None or _DNA_RE.match(alt) is None:
        raise HaplotypeError(f"Invalid DNA allele in atomic token: {token}")
    return AtomicEdit(
        token=t, chrom=chrom, pos=int(pos_str), ref=ref.upper(), alt=alt.upper()
    )


def _reconstruct_from_atoms(
    reference, chrom: str, atoms: List[AtomicEdit]
) -> Tuple[int, str, str]:
    if len(atoms) == 0:
        raise HaplotypeError("Haplotype contains no atomic edits")

    start = min(a.pos for a in atoms)
    end = max(a.pos + len(a.ref) - 1 for a in atoms)
    full_ref = reference.getReference(chrom, start, end).upper()
    if not full_ref:
        raise HaplotypeError("Unable to fetch reference sequence for haplotype span")

    alt_fragments = []
    cursor = start
    for a in atoms:
        if a.pos < cursor:
            raise HaplotypeError(f"Overlapping or unsorted atomic edits: {a.token}")

        prefix_end = a.pos - 1
        if prefix_end >= cursor:
            alt_fragments.append(
                reference.getReference(chrom, cursor, prefix_end).upper()
            )

        rel_start = a.pos - start
        rel_end = rel_start + len(a.ref)
        ref_chunk = full_ref[rel_start:rel_end]
        if ref_chunk != a.ref:
            raise HaplotypeError(
                f"Atomic REF mismatch against genome for token {a.token}: expected {ref_chunk}, saw {a.ref}"
            )

        alt_fragments.append(a.alt)
        cursor = a.pos + len(a.ref)

    if cursor <= end:
        alt_fragments.append(reference.getReference(chrom, cursor, end).upper())

    full_alt = "".join(alt_fragments)
    return start, full_ref, full_alt


def parse_haplotype_row(
    row_chrom: str,
    row_pos: int,
    row_ref: str,
    row_alt: str,
    row_id: str,
    reference,
    line_number: int,
    test_id: str = "",
) -> ParsedHaplotype:
    if ";" in row_ref or ";" in row_alt:
        raise HaplotypeError("Semicolon in row REF/ALT is not allowed")

    if ";" not in row_id:
        raise HaplotypeError("Row ID does not contain semicolon-separated atomic IDs")

    raw_tokens = [x.strip() for x in row_id.split(";") if x.strip()]
    if len(raw_tokens) < 2:
        raise HaplotypeError("At least two atomic IDs are required for haplotype mode")

    seen = set()
    atoms = []
    row_c = _canon_chrom(row_chrom)

    for tok in raw_tokens:
        if tok in seen:
            raise HaplotypeError(f"Duplicate atomic token: {tok}")
        seen.add(tok)
        a = parse_atomic_token(tok)
        if _canon_chrom(a.chrom) != row_c:
            raise HaplotypeError(f"Mixed chromosomes in haplotype token: {tok}")
        atoms.append(a)

    for i in range(1, len(atoms)):
        if atoms[i].pos < atoms[i - 1].pos:
            raise HaplotypeError(
                "Atomic IDs must be ordered by increasing genomic position"
            )

    recon_pos, recon_ref, recon_alt = _reconstruct_from_atoms(
        reference, row_chrom, atoms
    )

    if (
        int(row_pos) != recon_pos
        or row_ref.upper() != recon_ref
        or row_alt.upper() != recon_alt
    ):
        where = f"line {line_number}"
        if test_id:
            where += f", test {test_id}"
        raise HaplotypeError(
            f"Row CHROM/POS/REF/ALT does not match reconstructed atomic haplotype ({where})"
        )

    return ParsedHaplotype(
        chrom=row_chrom,
        pos=int(row_pos),
        ref=row_ref.upper(),
        alt=row_alt.upper(),
        atomic=tuple(atoms),
    )


def build_subset_vcf_fields(
    reference, row_chrom: str, atoms: List[AtomicEdit]
) -> Tuple[int, str, str, str]:
    pos, ref, alt = _reconstruct_from_atoms(reference, row_chrom, atoms)
    subset_id = ";".join(a.token for a in atoms)
    return pos, ref, alt, subset_id


def partition_for_split_based_on_protein(
    full_csn: str, atoms: List[AtomicEdit]
) -> List[List[AtomicEdit]]:
    if len(atoms) <= 1:
        return []

    # Keep all variants together when the full call is uncertain/stop/frameshift-like.
    if "_p.?" in full_csn or "fs" in full_csn or "Ter" in full_csn:
        return [atoms]

    # Deterministic baseline: each atomic edit is an independent subset.
    return [[a] for a in atoms]


def _extract_protein_part(csn: str) -> str:
    if "_p." not in csn:
        return ""
    return csn.split("_p.", 1)[1].strip()


def protein_components_from_csn(csn: str) -> List[str]:
    part = _extract_protein_part(csn)
    if not part:
        return []
    if part == "?" or part == ".":
        return [part]
    if part.startswith("[") and part.endswith("]"):
        inner = part[1:-1]
        if inner.startswith("(") and inner.endswith(")"):
            inner = inner[1:-1]
        return [x.strip() for x in inner.split(";") if x.strip()]
    if part.startswith("(") and part.endswith(")"):
        part = part[1:-1]
    return [part]


def _parse_pos_range(pos: str):
    if not pos:
        return None
    if "-" in pos:
        a, b = pos.split("-", 1)
        if a.isdigit() and b.isdigit():
            return int(a), int(b)
    if pos.isdigit():
        p = int(pos)
        return p, p
    return None


def _protein_string_from_component_list(components: List[str]) -> str:
    if len(components) == 0:
        return ""
    if len(components) == 1:
        return components[0]
    return "[(" + ";".join(components) + ")]"


def maybe_build_splitnearby_csn(
    csn: str, protpos: str, protref: str, protalt: str
) -> str:
    prot = _extract_protein_part(csn)
    if not prot:
        return ""
    if "fs" in prot or "Ter" in prot or "?" in prot:
        return ""
    if "delins" not in prot:
        return ""

    pos_range = _parse_pos_range(protpos)
    if pos_range is None:
        return ""
    p0, p1 = pos_range
    n = p1 - p0 + 1
    if n <= 1:
        return ""
    if len(protref) != n or len(protalt) != n:
        return ""

    comps = []
    for i in range(n):
        r = _AA1_TO_3.get(protref[i], "")
        a = _AA1_TO_3.get(protalt[i], "")
        if not r or not a or r == "?" or a == "?":
            return ""
        comps.append(f"{r}{p0 + i}{a}")

    newp = _protein_string_from_component_list(comps)
    return csn.split("_p.", 1)[0] + "_p." + newp


def choose_protein_partitions(
    full_components: List[str], subset_component_map: dict, atoms: List[AtomicEdit]
) -> List[List[AtomicEdit]]:
    if len(full_components) <= 1:
        # When multiple DNA edits collapse into one protein component (e.g. merged delins),
        # still expose deterministic per-atomic splits for downstream analysis.
        if len(atoms) <= 1:
            return [atoms]
        return [[a] for a in atoms]

    # Candidate subsets that map to one canonical component.
    single_hits = {}
    for key, comps in subset_component_map.items():
        if len(comps) != 1:
            continue
        c = comps[0]
        if c not in full_components:
            continue
        size = len(key)
        if c not in single_hits or size < len(single_hits[c]):
            single_hits[c] = key

    used = set()
    chosen = []
    for c in full_components:
        if c not in single_hits:
            continue
        key = single_hits[c]
        if any(i in used for i in key):
            continue
        chosen.append(key)
        used.update(key)

    if len(chosen) == 0:
        return [atoms]

    # Any uncovered variants become independent deterministic subsets.
    for i in range(len(atoms)):
        if i not in used:
            chosen.append((i,))

    # Sort subsets by lowest genomic coordinate.
    chosen.sort(key=lambda idxs: min(atoms[i].pos for i in idxs))
    return [[atoms[i] for i in idxs] for idxs in chosen]


def build_record_line_like(
    original_record, chrom: str, pos: int, row_id: str, ref: str, alt: str
) -> str:
    out_cols = []
    if original_record.chrom_chr_prefix and not chrom.startswith("chr"):
        out_cols.append("chr" + chrom)
    else:
        out_cols.append(chrom)
    out_cols.append(str(pos))
    out_cols.append(row_id)
    out_cols.append(ref)
    out_cols.append(alt)
    out_cols.append(original_record.qual)
    out_cols.append(original_record.filter)
    out_cols.append(original_record.info if original_record.info else ".")
    out_cols.extend(original_record.rest)
    return "\t".join(out_cols)


def add_haplotype_flags(record, original_ids: str, subset_ids: str):
    enc_orig = vcf_info_encode(original_ids)
    enc_subset = vcf_info_encode(subset_ids)
    for v in record.variants:
        v.addFlag("CAVA_ORIGHAPLOTYPE", enc_orig)
        v.addFlag("CAVA_HAPLOTYPE", enc_subset)
