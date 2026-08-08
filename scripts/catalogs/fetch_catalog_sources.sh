#!/usr/bin/env bash
set -euo pipefail

workspace="${1:-catalog-build}"
mkdir -p "$workspace"/{downloads,gtf,fasta,transcripts,reports}

validate_download() {
  local path="$1"
  [[ -s "$path" ]] || { echo "Downloaded file is empty: $path" >&2; return 1; }
  if [[ "$path" == *.gz ]]; then
    gzip -t "$path"
  fi
}

fetch() {
  local url="$1" output="$2"
  if [[ -s "$output" ]]; then
    if validate_download "$output"; then
      echo "Using existing $output"
      return
    fi
    echo "WARNING: existing file is invalid, re-downloading: $output" >&2
    rm -f "$output"
  fi
  if [[ -e "$output.part" ]] && ! validate_download "$output.part" >/dev/null 2>&1; then
    echo "WARNING: stale partial download is invalid, removing: $output.part" >&2
    rm -f "$output.part"
  fi
  echo "Downloading $url"
  if ! curl --fail --location --retry 4 --retry-delay 3 --continue-at - \
    --output "$output.part" "$url"; then
    echo "WARNING: curl download failed for $url; falling back to python urllib" >&2
    python3 - "$url" "$output.part" <<'PY_DL'
import shutil
import sys
import urllib.request

url, out = sys.argv[1], sys.argv[2]
with urllib.request.urlopen(url, timeout=180) as resp, open(out, "wb") as fh:
    shutil.copyfileobj(resp, fh)
PY_DL
  fi
  validate_download "$output.part"
  mv "$output.part" "$output"
}

# MANE 1.5 (current as of July 2026).
mane_base="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5"
fetch "$mane_base/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz" \
  "$workspace/gtf/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"
fetch "$mane_base/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz" \
  "$workspace/gtf/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
fetch "$mane_base/MANE.GRCh38.v1.5.refseq_rna.fna.gz" \
  "$workspace/downloads/MANE.GRCh38.v1.5.refseq_rna.fna.gz"
fetch "$mane_base/MANE.GRCh38.v1.5.ensembl_rna.fna.gz" \
  "$workspace/downloads/MANE.GRCh38.v1.5.ensembl_rna.fna.gz"
fetch "$mane_base/MANE.GRCh38.v1.5.summary.txt.gz" \
  "$workspace/downloads/MANE.GRCh38.v1.5.summary.txt.gz"
# NCBI GRCh38.p14 top-level assembly, including alternate loci and patches.
fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
  "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  "$workspace/fasta/GCF_000001405.40_GRCh38.p14_assembly_report.txt"

# Current Ensembl GRCh38, release 116.
ens116="https://ftp.ensembl.org/pub/release-116"
fetch "$ens116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz" \
  "$workspace/gtf/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz"
fetch "$ens116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" \
  "$workspace/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
fetch "$ens116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
  "$workspace/downloads/Homo_sapiens.GRCh38.116.cdna.all.fa.gz"

# Stable Ensembl GRCh37 annotation, release 75.
ens75="https://ftp.ensembl.org/pub/release-75"
fetch "$ens75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz" \
  "$workspace/gtf/Homo_sapiens.GRCh37.75.gtf.gz"
fetch "$ens75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz" \
  "$workspace/fasta/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz"
fetch "$ens75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz" \
  "$workspace/downloads/Homo_sapiens.GRCh37.75.cdna.all.fa.gz"

# Current GENCODE GRCh38, release 50. Use comprehensive annotation and the
# complete genome, not primary-only files, so alternate loci/patches including
# the sequence carrying GSTT1 are available.
gencode50="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50"
fetch "$gencode50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz" \
  "$workspace/gtf/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"
fetch "$gencode50/GRCh38.p14.genome.fa.gz" \
  "$workspace/fasta/GRCh38.p14.genome.fa.gz"
fetch "$gencode50/gencode.v50.transcripts.fa.gz" \
  "$workspace/downloads/gencode.v50.transcripts.fa.gz"
fetch "$gencode50/gencode.v50.metadata.RefSeq.gz" \
  "$workspace/downloads/gencode.v50.metadata.RefSeq.gz"

# NCBI's authoritative GRCh37.p13 RefSeq mapping. Current mapped annotation
# includes MANE tags and is the first choice for the experimental GRCh37 RefSeq
# catalog; it is safer than independently lifting GRCh38 exon coordinates.
gr37_ncbi="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13"
fetch "$gr37_ncbi/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" \
  "$workspace/gtf/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
fetch "$gr37_ncbi/GCF_000001405.25_GRCh37.p13_genomic.fna.gz" \
  "$workspace/fasta/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"
fetch "$gr37_ncbi/GCF_000001405.25_GRCh37.p13_rna.fna.gz" \
  "$workspace/downloads/GCF_000001405.25_GRCh37.p13_rna.fna.gz"
fetch "$gr37_ncbi/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
  "$workspace/fasta/GCF_000001405.25_GRCh37.p13_assembly_report.txt"

python3 "$(dirname "$0")/make_transcript_lists.py" --workspace "$workspace"

python3 - "$workspace" <<'PY_CHECKSUMS'
import hashlib
import sys
from pathlib import Path

workspace = Path(sys.argv[1]).resolve()
output = workspace / "reports/source_files.sha256"
files = sorted(
    path
    for directory in (workspace / "downloads", workspace / "gtf", workspace / "fasta")
    for path in directory.rglob("*")
    if path.is_file() and not path.name.endswith(".part")
)
with output.open("w", encoding="utf-8", newline="\n") as handle:
    for path in files:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                digest.update(block)
        handle.write("{}  {}\n".format(digest.hexdigest(), path.relative_to(workspace)))
print(output)
PY_CHECKSUMS

echo "Source files and transcript lists are ready under $workspace"
