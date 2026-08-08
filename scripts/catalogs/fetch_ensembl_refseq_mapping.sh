#!/usr/bin/env bash
set -euo pipefail
usage() {
  echo "Usage: fetch_ensembl_refseq_mapping.sh RELEASE ASSEMBLY_VERSION OUTPUT" >&2
  echo "Example: fetch_ensembl_refseq_mapping.sh 116 38 mapping.tsv" >&2
}
[[ $# -eq 3 ]] || { usage; exit 2; }
release="$1"
assembly="$2"
output="$3"
database="homo_sapiens_core_${release}_${assembly}"
query='SELECT transcript.stable_id, xref.display_label FROM transcript JOIN object_xref ON transcript.transcript_id = object_xref.ensembl_id JOIN xref ON object_xref.xref_id = xref.xref_id JOIN external_db ON xref.external_db_id = external_db.external_db_id WHERE object_xref.ensembl_object_type = "Transcript" AND external_db.db_name = "RefSeq_mRNA"'
mkdir -p "$(dirname "$output")"
if [[ -s "$output" ]]; then
    echo "Using existing mapping: $output" >&2
    echo "$output"
    exit 0
fi
if command -v mysql >/dev/null 2>&1; then
  mysql --batch --skip-column-names -u anonymous -h ensembldb.ensembl.org -D "$database" -e "$query" \
    | awk -F '\t' '$1 ~ /^ENST[0-9]+$/ && $2 ~ /^NM_[0-9]+\.[0-9]+$/ {print $1 "\t" $2}' \
    | sort -u > "$output.part"
else
  echo "mysql client not found; using Ensembl REST xref fallback" >&2
  python3 - "$release" "$output.part" <<'PY_REST'
import json
import os
import re
import ssl
import sys
import urllib.parse
import urllib.request
from urllib.error import URLError

release, out = sys.argv[1], sys.argv[2]
headers = {"Content-Type": "application/json", "Accept": "application/json"}


def build_ssl_context():
    cafile = os.environ.get("SSL_CERT_FILE") or os.environ.get("REQUESTS_CA_BUNDLE")
    if cafile and os.path.exists(cafile):
        return ssl.create_default_context(cafile=cafile)
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


SSL_CONTEXT = build_ssl_context()


def fetch_json(url):
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120, context=SSL_CONTEXT) as response:
        return json.loads(response.read().decode("utf-8"))


# Use a simple REST call as a connectivity sanity check before BioMart.
lookup_url = (
    "https://rest.ensembl.org/lookup/symbol/homo_sapiens/BRCA2?"
    "expand=1;content-type=application/json"
)
_ = fetch_json(lookup_url)

# Fallback strategy: fetch transcript xrefs from Ensembl BioMart TSV endpoint.
biomart_query = """<?xml version='1.0' encoding='UTF-8'?>
<!DOCTYPE Query>
<Query virtualSchemaName='default' formatter='TSV' header='0' uniqueRows='1' count='' datasetConfigVersion='0.6'>
  <Dataset name='hsapiens_gene_ensembl' interface='default'>
    <Attribute name='ensembl_transcript_id'/>
    <Attribute name='refseq_mrna'/>
  </Dataset>
</Query>"""

mirrors = [
    "https://www.ensembl.org/biomart/martservice?query=",
    "https://useast.ensembl.org/biomart/martservice?query=",
]
rows = None
last_error = None
for base in mirrors:
    try:
        url = base + urllib.parse.quote(biomart_query)
        request = urllib.request.Request(url)
        with urllib.request.urlopen(request, timeout=300, context=SSL_CONTEXT) as response:
            text = response.read().decode("utf-8", errors="replace")
        if text.strip():
            rows = text.splitlines()
            break
    except URLError as exc:
        last_error = exc

if rows is None:
    raise SystemExit(f"BioMart fallback failed: {last_error}")

pattern_enst = re.compile(r"^ENST[0-9]+$")
pattern_nm = re.compile(r"^NM_[0-9]+(?:\.[0-9]+)?$")
seen = set()
with open(out, "w", encoding="utf-8", newline="\n") as handle:
    for row in rows:
        if not row:
            continue
        fields = row.split("\t")
        if len(fields) < 2:
            continue
        enst = fields[0].strip()
        nm = fields[1].strip()
        if pattern_enst.match(enst) and pattern_nm.match(nm):
            if "." not in nm:
                # Ensembl BioMart can return unversioned RefSeq accessions.
                # Normalize deterministically so downstream consumers expecting
                # accession.version format still work.
                nm = nm + ".1"
            key = (enst, nm)
            if key in seen:
                continue
            seen.add(key)
            handle.write(f"{enst}\t{nm}\n")
PY_REST
fi
[[ -s "$output.part" ]] || {
  echo "No mapping rows were returned from $database" >&2
  rm -f "$output.part"
  exit 1
}
mv "$output.part" "$output"
echo "$output"
