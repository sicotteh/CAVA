#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel 2>/dev/null)" || {
  echo "Run this script inside the checked-out CAVA Git repository." >&2
  exit 2
}
cd "$repo_root"

for required in cava/ensembldb/README.md config_template.txt pyproject.toml cava_catalogs.tsv; do
  [[ -f "$required" ]] || {
    echo "Expected CAVA file is missing: $required" >&2
    exit 2
  }
done

git lfs version >/dev/null 2>&1 || {
  echo "git-lfs is required. Install it, run 'git lfs install', and rerun." >&2
  exit 2
}

git lfs install
mkdir -p catalogs/files
touch catalogs/files/.gitkeep

# Catalog payloads are immutable data objects. The manifest remains ordinary
# Git text so users can list catalogs without downloading any LFS content.
git lfs track "catalogs/files/*.gz"
git lfs track "catalogs/files/*.gz.tbi"
git lfs track "catalogs/files/*.txt"
git lfs track "catalogs/files/*.cesis"

git add .gitattributes cava_catalogs.tsv catalogs/files/.gitkeep

git check-attr filter diff merge text -- \
  catalogs/files/example.gz \
  catalogs/files/example.gz.tbi \
  catalogs/files/example.txt \
  catalogs/files/example.cesis

cat <<'EOF'
Git LFS is configured for CAVA catalog payloads.

Review and commit the one-time setup:
  git diff --cached -- .gitattributes cava_catalogs.tsv
  git commit -m "Configure Git LFS for CAVA catalogs"
  git push origin HEAD

Code-only commits do not require rebuilding, staging, or re-uploading catalogs.
EOF
