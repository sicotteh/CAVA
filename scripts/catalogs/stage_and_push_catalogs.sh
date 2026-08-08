#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: stage_and_push_catalogs.sh [-m MESSAGE] [--no-push] CATALOG_NAME [...]

Stages only the selected catalog payloads and cava_catalogs.tsv. It does not
create a GitHub Release and is never needed for an unrelated code-only commit.
EOF
}
message="Update CAVA transcript catalogs"
push=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--message) message="$2"; shift 2 ;;
    --no-push) push=0; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    *) break ;;
  esac
done
[[ $# -ge 1 ]] || { usage >&2; exit 2; }
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"
git lfs version >/dev/null
python3 scripts/catalogs/validate_manifest.py --repo "$repo_root" --check-files

paths=(cava_catalogs.tsv)
for name in "$@"; do
  while IFS= read -r path; do
    [[ -n "$path" ]] && paths+=("$path")
  done < <(python3 - "$name" <<'PY'
import csv
import sys
with open('cava_catalogs.tsv', newline='', encoding='utf-8') as handle:
    rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['catalog_name'] == sys.argv[1]]
if len(rows) != 1:
    raise SystemExit("Expected one manifest row for {!r}; found {}".format(sys.argv[1], len(rows)))
for key in ('lfs_path', 'index_path', 'transcript_map_path', 'selenocysteine_path'):
    print(rows[0][key])
PY
)
done

git add -- "${paths[@]}"
python3 - "$@" <<'PY'
import csv
import re
import subprocess
import sys

names = sys.argv[1:]
with open('cava_catalogs.tsv', newline='', encoding='utf-8') as handle:
    rows = {row['catalog_name']: row for row in csv.DictReader(handle, delimiter='\t')}

groups = (
    ('lfs_path', 'lfs_oid_sha256', 'lfs_size'),
    ('index_path', 'index_lfs_oid_sha256', 'index_lfs_size'),
    ('transcript_map_path', 'transcript_map_lfs_oid_sha256', 'transcript_map_lfs_size'),
    ('selenocysteine_path', 'selenocysteine_lfs_oid_sha256', 'selenocysteine_lfs_size'),
)
for name in names:
    row = rows[name]
    for path_key, oid_key, size_key in groups:
        path = row[path_key]
        pointer = subprocess.check_output(['git', 'show', ':' + path], text=True)
        fields = dict(line.split(' ', 1) for line in pointer.strip().splitlines() if ' ' in line)
        if fields.get('version') != 'https://git-lfs.github.com/spec/v1':
            raise SystemExit('Staged object is not a Git LFS pointer: ' + path)
        match = re.fullmatch(r'sha256:([0-9a-f]{64})', fields.get('oid', ''))
        if not match or match.group(1) != row[oid_key]:
            raise SystemExit('Staged pointer OID mismatch for ' + path)
        if fields.get('size') != row[size_key]:
            raise SystemExit('Staged pointer size mismatch for ' + path)
PY

git diff --cached --stat
git commit -m "$message" -- "${paths[@]}"
if [[ "$push" -eq 1 ]]; then
  # The LFS pre-push hook transfers only objects that are not already remote.
  git push origin HEAD
else
  echo "Committed locally; --no-push was requested."
fi
