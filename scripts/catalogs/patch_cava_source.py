#!/usr/bin/env python3
"""Apply the small source changes required by the catalog workflow.

The edits target the current CAVA 2.0.15 master layout and are idempotent:

* add the ``cava_data`` console script and install the two root data files;
* teach ``cava/MANE.py`` to accept explicitly selected alt-contig genes or
  transcript accessions;
* pass that selection into ``mane_db_prep.py`` so GSTT1 can be retained while
  the default primary-chromosome behavior remains unchanged;
* make CAVA use a same-basename ``.cesis`` sidecar when one is installed next
  to the selected transcript catalog.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


class PatchError(RuntimeError):
    pass


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise PatchError(
            "Expected exactly one {} anchor, found {}. The source may have changed.".format(
                label, count
            )
        )
    return text.replace(old, new, 1)


def write_if_changed(path: Path, original: str, updated: str) -> bool:
    if original == updated:
        return False
    backup = path.with_name(path.name + ".catalog-backup")
    if not backup.exists():
        backup.write_text(original, encoding="utf-8")
    path.write_text(updated, encoding="utf-8")
    return True


def patch_pyproject(path: Path) -> bool:
    original = path.read_text(encoding="utf-8")
    updated = original

    script_line = 'cava_data = "cava.cava_data:main"'
    if script_line not in updated:
        if "[project.scripts]" in updated:
            updated = replace_once(
                updated,
                "[project.scripts]\n",
                "[project.scripts]\n{}\n".format(script_line),
                "[project.scripts]",
            )
        else:
            updated = updated.rstrip() + "\n\n[project.scripts]\n{}\n".format(script_line)

    data_line = '"share/cava" = ["config_template.txt", "cava_catalogs.tsv"]'
    if data_line not in updated:
        if "[tool.setuptools.data-files]" in updated:
            updated = replace_once(
                updated,
                "[tool.setuptools.data-files]\n",
                "[tool.setuptools.data-files]\n{}\n".format(data_line),
                "[tool.setuptools.data-files]",
            )
        else:
            updated = (
                updated.rstrip()
                + "\n\n[tool.setuptools.data-files]\n{}\n".format(data_line)
            )
    return write_if_changed(path, original, updated)


def patch_mane_cli(path: Path) -> bool:
    original = path.read_text(encoding="utf-8")
    if "--include-alt-gene" in original and "_split_csv_options" in original:
        return False
    updated = original

    import_anchor = "import os\n"
    helper = '''import os\n\n\ndef _split_csv_options(values):\n    result = set()\n    for value in values or []:\n        for item in value.split(','):\n            item = item.strip()\n            if item:\n                result.add(item)\n    return result\n'''
    updated = replace_once(updated, import_anchor, helper, "MANE.py import")

    option_anchor = '''parser.add_option("--no_hg19",  action='store_false', default=True, dest='no_hg19', help="Set this to skip hg19 builds")\n'''
    option_replacement = option_anchor + '''parser.add_option(\n    "--include-alt-gene",\n    action="append",\n    default=[],\n    dest="include_alt_genes",\n    help=(\n        "Gene symbol to retain when its MANE annotation is on an alternate, "\n        "fix, patch or scaffold contig. Repeat the option or use commas. "\n        "Example: --include-alt-gene GSTT1"\n    ),\n)\nparser.add_option(\n    "--include-alt-transcript",\n    action="append",\n    default=[],\n    dest="include_alt_transcripts",\n    help=(\n        "Transcript accession to retain on a non-primary contig. Repeat the "\n        "option or use commas. A versioned accession also matches MANE's "\n        "location suffix form, for example NM_000853.4_1."\n    ),\n)\n'''
    updated = replace_once(
        updated, option_anchor, option_replacement, "MANE.py alt-contig options"
    )

    parse_anchor = "(options, args) = parser.parse_args()\n\noptions.select = False\n"
    parse_replacement = '''(options, args) = parser.parse_args()\n\noptions.include_alt_genes = _split_csv_options(options.include_alt_genes)\noptions.include_alt_transcripts = _split_csv_options(options.include_alt_transcripts)\noptions.select = False\n'''
    updated = replace_once(
        updated, parse_anchor, parse_replacement, "MANE.py parsed options"
    )
    return write_if_changed(path, original, updated)


def patch_mane_builder(path: Path) -> bool:
    original = path.read_text(encoding="utf-8")
    if "def _mane_record_is_selected" in original and "options=options" in original:
        return False
    updated = original

    if "import re\n" not in updated:
        updated = replace_once(updated, "import os\n", "import os\nimport re\n", "mane_db_prep import")

    helper_anchor = "failed_conversions['ENST'] = set()\n\n\n"
    helper = r'''failed_conversions['ENST'] = set()

PRIMARY_MANE_CONTIGS = {
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
    'MT', 'X', 'Y'
}
reported_alt_mane_records = set()


def _mane_transcript_base(identifier):
    """Remove only MANE's trailing alternate-location suffix.

    RefSeq version numbers and Ensembl version numbers remain part of the
    accession. For example ``NM_000853.4_1`` becomes ``NM_000853.4``.
    """
    return re.sub(r'_[0-9]+$', '', identifier or '')


def _mane_record_is_selected(contig, tags, options):
    if contig in PRIMARY_MANE_CONTIGS:
        return True

    transcript = getValue(tags, 'transcript_id') or ''
    gene_values = {
        value.upper()
        for value in (
            getValue(tags, 'gene_name'),
            getValue(tags, 'gene'),
            getValue(tags, 'gene_id'),
        )
        if value
    }
    requested_genes = {
        value.upper() for value in getattr(options, 'include_alt_genes', set())
    }
    requested_transcripts = set(
        getattr(options, 'include_alt_transcripts', set())
    )
    selected = bool(gene_values & requested_genes)
    selected = selected or transcript in requested_transcripts
    selected = selected or _mane_transcript_base(transcript) in requested_transcripts
    if selected:
        key = (contig, transcript, tuple(sorted(gene_values)))
        if key not in reported_alt_mane_records:
            reported_alt_mane_records.add(key)
            print(
                'Including requested MANE alt-contig transcript: '
                + contig + '\t' + ','.join(sorted(gene_values)) + '\t' + transcript
            )
    return selected


'''
    updated = replace_once(
        updated, helper_anchor, helper, "mane_db_prep alt-contig helper"
    )

    updated = replace_once(
        updated,
        "def parse_GTF(filename='', genesdata=None):\n",
        "def parse_GTF(filename='', genesdata=None, options=None):\n",
        "mane_db_prep parse_GTF signature",
    )

    start_marker = "        cols = line.split('\\t')\n\n        # Only consider transcripts on the following chromosomes\n"
    start = updated.find(start_marker)
    if start < 0:
        raise PatchError("Could not find MANE chromosome-filter block")
    end_marker = "        tags = cols[8].split(';')\n"
    end = updated.find(end_marker, start)
    if end < 0:
        raise PatchError("Could not find MANE tag-parsing block")
    end += len(end_marker)
    replacement = '''        cols = line.split('\\t')\n\n        # Consider only records used to build a CAVA transcript.\n        if cols[2] not in ['exon', 'transcript', 'start_codon', 'stop_codon']:\n            continue\n\n        # Parse attributes before filtering the chromosome. This preserves the\n        # historical primary-contig default while allowing a maintainer to\n        # explicitly retain selected alternate-locus transcripts such as GSTT1.\n        tags = cols[8].split(';')\n        if not _mane_record_is_selected(cols[0], tags, options):\n            continue\n'''
    updated = updated[:start] + replacement + updated[end:]

    call_pattern = re.compile(
        r"    transcript, prevenst, first, genesdata = parse_GTF\("
        r"filename=source_compressed_gtf,\n\s+genesdata=genesdata\)\n"
    )
    call_new = '''    transcript, prevenst, first, genesdata = parse_GTF(\n        filename=source_compressed_gtf,\n        genesdata=genesdata,\n        options=options,\n    )\n'''
    updated, count = call_pattern.subn(call_new, updated, count=1)
    if count != 1:
        raise PatchError("Could not find mane_db_prep parse_GTF call")
    return write_if_changed(path, original, updated)


def patch_selenodata_discovery(path: Path) -> bool:
    original = path.read_text(encoding="utf-8")
    if "adjacent_selenofile" in original:
        return False
    old = '''        else:\n            fid = _open_ensembldb_resource("SECIS_in_refseq_pos.txt")\n'''
    new = '''        else:\n            # A cava_data installation places the matching SECIS/selenocysteine\n            # annotation next to the catalog. This keeps config generation\n            # backward compatible: only @ensembl needs to change.\n            catalog_file = str(options.args.get('ensembl', ''))\n            adjacent_selenofile = re.sub(r'\\.gz$', '.cesis', catalog_file)\n            if adjacent_selenofile != catalog_file and os.path.isfile(adjacent_selenofile):\n                fid = open(adjacent_selenofile, 'r')\n            else:\n                fid = _open_ensembldb_resource("SECIS_in_refseq_pos.txt")\n'''
    updated = replace_once(original, old, new, "same-basename .cesis discovery")
    return write_if_changed(path, original, updated)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path, nargs="?", default=Path.cwd())
    args = parser.parse_args()
    repo = args.repo.resolve()

    required = [
        repo / "pyproject.toml",
        repo / "cava" / "MANE.py",
        repo / "cava" / "ensembldb" / "mane_db_prep.py",
        repo / "cava" / "utils" / "data.py",
    ]
    for path in required:
        if not path.is_file():
            raise PatchError("Missing expected CAVA file: {}".format(path))

    changed = []
    operations = (
        (required[0], patch_pyproject),
        (required[1], patch_mane_cli),
        (required[2], patch_mane_builder),
        (required[3], patch_selenodata_discovery),
    )
    for path, operation in operations:
        if operation(path):
            changed.append(path.relative_to(repo))

    if changed:
        print("Patched:")
        for path in changed:
            print("  {}".format(path))
    else:
        print("CAVA catalog source changes are already present.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PatchError as exc:
        raise SystemExit("patch_cava_source.py: error: {}".format(exc))
