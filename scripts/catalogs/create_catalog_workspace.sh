#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel 2>/dev/null)" || {
  echo "Run inside the CAVA repository." >&2
  exit 2
}
workspace="${1:-$repo_root/catalog-build}"
mkdir -p "$workspace"/{downloads,gtf,fasta,transcripts,output,logs,reports,tmp}
printf '%s\n' "$workspace"
