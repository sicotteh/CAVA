#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: build_catalogs.sh [--include-experimental] [WORKSPACE]

Builds CAVA catalogs from the current checked-out master source. The stable
MANE catalogs are built through cava/MANE.py; the source patch adds explicit
alt-contig selection and this script requests GSTT1.

The script updates local catalog files and cava_catalogs.tsv but does not stage,
commit, push, tag, or create a GitHub Release.
EOF
}
include_experimental=0
if [[ "${1:-}" == "--include-experimental" ]]; then
  include_experimental=1
  shift
fi
[[ "${1:-}" != "-h" && "${1:-}" != "--help" ]] || { usage; exit 0; }

repo_root="$(git rev-parse --show-toplevel)"
workspace="${1:-$repo_root/catalog-build}"
script_dir="$repo_root/scripts/catalogs"
mkdir -p "$workspace"/{build,downloads,gtf,fasta,prepared-fasta,transcripts,validated-transcripts,validated-gtf,mappings,output,logs,reports,tmp,selenodata}

for required in \
  cava/MANE.py \
  cava/EnsemblDB.py \
  cava/RefSeqDB.py \
  cava/ensembldb/create_selenodata.py \
  cava/ensembldb/README.md \
  config_template.txt \
  cava_catalogs.tsv; do
  [[ -f "$repo_root/$required" ]] || { echo "Missing $repo_root/$required" >&2; exit 2; }
done

python3 "$repo_root/cava/MANE.py" -h 2>&1 | grep -q -- '--include-alt-gene' || {
  echo "CAVA MANE.py has not been patched for explicit alt-contig inclusion." >&2
  echo "Run: python3 scripts/catalogs/patch_cava_source.py '$repo_root'" >&2
  exit 2
}

"$script_dir/fetch_catalog_sources.sh" "$workspace"

prepared() {
  local source="$1"
  local base
  base="$(basename "$source")"
  base="${base%.gz}"
  printf '%s/prepared-fasta/%s' "$workspace" "$base"
}

ensure_prepared() {
  local source="$1"
  python3 "$script_dir/prepare_reference_fastas.py" --workspace "$workspace" --source "$source" >&2
}

drop_prepared() {
  local source="$1"
  local prepared_path
  prepared_path="$(prepared "$source")"
  rm -f "$prepared_path" "$prepared_path.fai" "$prepared_path.aliases.tsv"
}

clean_validated_ids() {
  local source="$1" target="$2"
  python3 - "$source" "$target" <<'PY'
import sys
from pathlib import Path
excluded = {"NM_001424184.1", "ENST00000434431.2"}
source, target = map(Path, sys.argv[1:])
values = sorted({line.strip() for line in source.read_text().splitlines() if line.strip()} - excluded)
target.write_text("".join(value + "\n" for value in values), encoding="utf-8")
print("{}: {} validated transcripts after explicit partial-transcript exclusions".format(target, len(values)))
PY
}

validate_gtf() {
  local name="$1" gtf="$2" genome="$3" rna="$4" list="$5" require_gstt1="$6" allow_mismatches="$7"
  local raw="$workspace/validated-transcripts/${name}.raw.txt"
  local final="$workspace/validated-transcripts/${name}.txt"
  local filtered_gtf="$workspace/validated-gtf/${name}.gtf.gz"
  local command=(
    python3 "$script_dir/validate_gtf_sequences.py"
    --gtf "$gtf"
    --genome-fasta "$genome"
    --rna-fasta "$rna"
    --transcripts "$list"
    --report "$workspace/reports/${name}.source-sequences.tsv"
    --accepted-transcripts "$raw"
    --filtered-gtf "$filtered_gtf"
  )
  [[ "$require_gstt1" == 1 ]] && command+=(--require-gene GSTT1)
  [[ "$allow_mismatches" == 1 ]] && command+=(--allow-mismatches)
  "${command[@]}" >&2
  clean_validated_ids "$raw" "$final" >&2
  printf '%s' "$final"
}

copy_products() {
  local built_catalog="$1" built_map="$2" final_base="$3"
  cp -f "$built_catalog" "$final_base.gz"
  cp -f "$built_catalog.tbi" "$final_base.gz.tbi"
  cp -f "$built_map" "$final_base.txt"
}

create_seleno_refseq() {
  local ids="$1" tag="$2" output="$3"
  run_seleno_with_retry \
    --repo "$repo_root" --refseqs "$ids" --tag "$tag" --output "$output"
}

create_seleno_ensembl() {
  local mapping="$1" tag="$2" output="$3"
  run_seleno_with_retry \
    --repo "$repo_root" --ensembl-refseq "$mapping" --tag "$tag" --output "$output"
}

enforce_protein_ids() {
  local name="$1" final_base="$2"
  python3 "$script_dir/enforce_companion_proteins.py" \
    --catalog "$final_base.gz" \
    --transcript-map "$final_base.txt" \
    --report "$workspace/reports/${name}.companion-proteins.tsv"
}

run_seleno_with_retry() {
  local max_attempts=3
  local attempt
  for ((attempt=1; attempt<=max_attempts; attempt++)); do
    if python3 "$script_dir/create_selenocysteine_sidecar.py" "$@"; then
      return 0
    fi
    echo "Selenocysteine sidecar generation failed (attempt $attempt/$max_attempts); retrying..." >&2
  done
  echo "Selenocysteine sidecar generation failed after $max_attempts attempts." >&2
  return 1
}

validate_final() {
  local name="$1" final_base="$2" genome="$3" rna="$4" require_gstt1="$5"
  local command=(
    python3 "$script_dir/validate_catalog_sequences.py"
    --catalog "$final_base.gz"
    --genome-fasta "$genome"
    --rna-fasta "$rna"
    --report "$workspace/reports/${name}.final-catalog-sequences.tsv"
  )
  [[ "$require_gstt1" == 1 ]] && command+=(--require-gene GSTT1)
  "${command[@]}"
  python3 "$script_dir/validate_catalog_companions.py" \
    --catalog "$final_base.gz" \
    --transcript-map "$final_base.txt" \
    --report "$workspace/reports/${name}.companion-map.tsv"
  gzip -t "$final_base.gz"
  [[ -s "$final_base.gz.tbi" && -s "$final_base.txt" && -s "$final_base.cesis" ]]
}

register_catalog() {
  local name="$1" final_base="$2" status="$3"
  python3 "$script_dir/update_catalog_manifest.py" "$name" \
    --repo "$repo_root" \
    --catalog "$final_base.gz" \
    --index "$final_base.gz.tbi" \
    --transcript-map "$final_base.txt" \
    --selenocysteine "$final_base.cesis" \
    --status "$status"
}

run_ensembl_builder() {
  local release="$1" build="$2" url="$3" source_gtf="$4" accepted="$5" output_name="$6" output_dir="$7"
  local run_dir="$workspace/tmp/${output_name}.run"
  mkdir -p "$run_dir" "$output_dir"
  local basename
  basename="$(basename "$url")"
  ln -sfn "$source_gtf" "$run_dir/$basename"
  (
    cd "$run_dir"
    PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}" \
      python3 "$repo_root/cava/EnsemblDB.py" \
      -e "$release" -b "$build" --no_hg19 \
      -u "$url" -i "$accepted" -o "$output_name" -D "$output_dir"
  ) 2>&1 | tee "$workspace/logs/${output_name}.log"
}

# ---------------------------------------------------------------------------
# MANE 1.5 GRCh38. Use the repository's MANE.py for both ID systems.

ensure_prepared "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
ensure_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"
mane_refseq_genome="$(prepared "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz")"
mane_ensembl_genome="$(prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz")"
mane_refseq_ids="$(validate_gtf \
  mane-1.5-grch38-refseq \
  "$workspace/gtf/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz" \
  "$mane_refseq_genome" \
  "$workspace/downloads/MANE.GRCh38.v1.5.refseq_rna.fna.gz" \
  "$workspace/transcripts/MANE_1.5_GRCh38_REFSEQ.txt" 1 1)"
mane_ensembl_ids="$(validate_gtf \
  mane-1.5-grch38-ensembl \
  "$workspace/gtf/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz" \
  "$mane_ensembl_genome" \
  "$workspace/downloads/MANE.GRCh38.v1.5.ensembl_rna.fna.gz" \
  "$workspace/transcripts/MANE_1.5_GRCh38_ENSEMBL.txt" 1 1)"

mane_out="$workspace/build/mane-1.5"
mane_run="$workspace/tmp/mane-1.5.run"
mkdir -p "$mane_out" "$mane_run"
cp -f "$workspace/validated-gtf/mane-1.5-grch38-refseq.gtf.gz" \
  "$mane_out/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"
cp -f "$workspace/validated-gtf/mane-1.5-grch38-ensembl.gtf.gz" \
  "$mane_out/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
(
  cd "$mane_run"
  PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}" \
    python3 "$repo_root/cava/MANE.py" \
    -e 1.5 -D "$mane_out" --no_hg19 --include-alt-gene GSTT1
) 2>&1 | tee "$workspace/logs/MANE_1.5.log"

mane_refseq_base="$workspace/output/CAVA_MANE_1.5_GRCh38_REFSEQ"
python3 "$script_dir/filter_catalog.py" \
  --catalog "$mane_out/MANE.GRCh38.v1.5.refseq_genomic.db.gz" \
  --transcript-map "$mane_out/MANE.GRCh38.v1.5.refseq_genomic.txt" \
  --allowed "$mane_refseq_ids" \
  --output-catalog "$mane_refseq_base.gz" \
  --output-map "$mane_refseq_base.txt" \
  --require-gene GSTT1
create_seleno_refseq "$mane_refseq_ids" mane_1_5_grch38_refseq "$mane_refseq_base.cesis"
enforce_protein_ids mane-1.5-grch38-refseq "$mane_refseq_base"
validate_final mane-1.5-grch38-refseq "$mane_refseq_base" "$mane_refseq_genome" \
  "$workspace/downloads/MANE.GRCh38.v1.5.refseq_rna.fna.gz" 1
register_catalog mane-1.5-grch38-refseq "$mane_refseq_base" published
if [[ "$include_experimental" != 1 ]]; then
  drop_prepared "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
fi

mane_mapping="$workspace/mappings/MANE_1.5.ensembl_refseq.tsv"
python3 "$script_dir/make_identifier_mappings.py" --mode mane \
  --input "$workspace/downloads/MANE.GRCh38.v1.5.summary.txt.gz" \
  --transcripts "$mane_ensembl_ids" --output "$mane_mapping"
mane_ensembl_base="$workspace/output/CAVA_MANE_1.5_GRCh38_ENSEMBL"
python3 "$script_dir/filter_catalog.py" \
  --catalog "$mane_out/MANE.GRCh38.v1.5.ensembl_genomic.db.gz" \
  --transcript-map "$mane_out/MANE.GRCh38.v1.5.ensembl_genomic.txt" \
  --allowed "$mane_ensembl_ids" \
  --output-catalog "$mane_ensembl_base.gz" \
  --output-map "$mane_ensembl_base.txt" \
  --require-gene GSTT1
create_seleno_ensembl "$mane_mapping" mane_1_5_grch38_ensembl "$mane_ensembl_base.cesis"
enforce_protein_ids mane-1.5-grch38-ensembl "$mane_ensembl_base"
validate_final mane-1.5-grch38-ensembl "$mane_ensembl_base" "$mane_ensembl_genome" \
  "$workspace/downloads/MANE.GRCh38.v1.5.ensembl_rna.fna.gz" 1
register_catalog mane-1.5-grch38-ensembl "$mane_ensembl_base" published
if [[ "$include_experimental" != 1 ]]; then
  drop_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"
fi

# ---------------------------------------------------------------------------
# Ensembl 116 GRCh38.

ensure_prepared "$workspace/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
ens116_genome="$(prepared "$workspace/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz")"
ens116_ids="$(validate_gtf \
  ensembl-116-grch38 \
  "$workspace/gtf/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz" \
  "$ens116_genome" \
  "$workspace/downloads/Homo_sapiens.GRCh38.116.cdna.all.fa.gz" \
  "$workspace/transcripts/Ensembl_116_GRCh38.txt" 1 1)"
ens116_url="https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz"
ens116_build="$workspace/build/ensembl-116"
run_ensembl_builder 116 GRCh38 "$ens116_url" \
  "$workspace/validated-gtf/ensembl-116-grch38.gtf.gz" \
  "$ens116_ids" CAVA_Ensembl_116_GRCh38 "$ens116_build"
ens116_base="$workspace/output/CAVA_Ensembl_116_GRCh38"
copy_products "$ens116_build/CAVA_Ensembl_116_GRCh38.gz" \
  "$ens116_build/CAVA_Ensembl_116_GRCh38.txt" "$ens116_base"
ens116_map="$workspace/mappings/ensembl_refseq_GRCh38_116.tsv"
"$script_dir/fetch_ensembl_refseq_mapping.sh" 116 38 "$ens116_map"
create_seleno_ensembl "$ens116_map" ensembl_grch38_116 "$ens116_base.cesis"
validate_final ensembl-116-grch38 "$ens116_base" "$ens116_genome" \
  "$workspace/downloads/Homo_sapiens.GRCh38.116.cdna.all.fa.gz" 1
register_catalog ensembl-116-grch38 "$ens116_base" published
drop_prepared "$workspace/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

# ---------------------------------------------------------------------------
# Native Ensembl 75 GRCh37.

ensure_prepared "$workspace/fasta/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz"
ens75_genome="$(prepared "$workspace/fasta/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz")"
ens75_ids="$(validate_gtf \
  ensembl-75-grch37 \
  "$workspace/gtf/Homo_sapiens.GRCh37.75.gtf.gz" \
  "$ens75_genome" \
  "$workspace/downloads/Homo_sapiens.GRCh37.75.cdna.all.fa.gz" \
  "$workspace/transcripts/Ensembl_75_GRCh37.txt" 0 1)"
ens75_url="https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
ens75_build="$workspace/build/ensembl-75"
run_ensembl_builder 75 GRCh37 "$ens75_url" \
  "$workspace/validated-gtf/ensembl-75-grch37.gtf.gz" \
  "$ens75_ids" CAVA_Ensembl_75_GRCh37 "$ens75_build"
ens75_base="$workspace/output/CAVA_Ensembl_75_GRCh37"
copy_products "$ens75_build/CAVA_Ensembl_75_GRCh37.gz" \
  "$ens75_build/CAVA_Ensembl_75_GRCh37.txt" "$ens75_base"
ens75_map="$workspace/mappings/ensembl_refseq_GRCh37_75.tsv"
"$script_dir/fetch_ensembl_refseq_mapping.sh" 75 37 "$ens75_map"
create_seleno_ensembl "$ens75_map" ensembl_grch37_75 "$ens75_base.cesis"
validate_final ensembl-75-grch37 "$ens75_base" "$ens75_genome" \
  "$workspace/downloads/Homo_sapiens.GRCh37.75.cdna.all.fa.gz" 0
register_catalog ensembl-75-grch37 "$ens75_base" published
if [[ "$include_experimental" != 1 ]]; then
  drop_prepared "$workspace/fasta/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz"
fi

# ---------------------------------------------------------------------------
# GENCODE 50 GRCh38 comprehensive ALL annotation.

ensure_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"
gencode_genome="$(prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz")"
gencode_ids="$(validate_gtf \
  gencode-50-grch38 \
  "$workspace/gtf/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz" \
  "$gencode_genome" \
  "$workspace/downloads/gencode.v50.transcripts.fa.gz" \
  "$workspace/transcripts/GENCODE_50_GRCh38.txt" 1 1)"
gencode_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"
gencode_build="$workspace/build/gencode-50"
run_ensembl_builder 116 GRCh38 "$gencode_url" \
  "$workspace/validated-gtf/gencode-50-grch38.gtf.gz" \
  "$gencode_ids" CAVA_GENCODE_50_GRCh38 "$gencode_build"
gencode_base="$workspace/output/CAVA_GENCODE_50_GRCh38"
copy_products "$gencode_build/CAVA_GENCODE_50_GRCh38.gz" \
  "$gencode_build/CAVA_GENCODE_50_GRCh38.txt" "$gencode_base"
gencode_map="$workspace/mappings/gencode_50_refseq.tsv"
python3 "$script_dir/make_identifier_mappings.py" --mode gencode \
  --input "$workspace/downloads/gencode.v50.metadata.RefSeq.gz" \
  --transcripts "$gencode_ids" --output "$gencode_map"
create_seleno_ensembl "$gencode_map" gencode_grch38_50 "$gencode_base.cesis"
validate_final gencode-50-grch38 "$gencode_base" "$gencode_genome" \
  "$workspace/downloads/gencode.v50.transcripts.fa.gz" 1
register_catalog gencode-50-grch38 "$gencode_base" published
drop_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"

# ---------------------------------------------------------------------------
# Experimental MANE 1.5 projection to GRCh37.

if [[ "$include_experimental" == 1 ]]; then
  ensure_prepared "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
  ensure_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"
  python3 "$script_dir/build_mane_grch37_candidates.py" \
    --mane-summary "$workspace/downloads/MANE.GRCh38.v1.5.summary.txt.gz" \
    --mane-refseq-grch38-gtf "$workspace/gtf/MANE.GRCh38.v1.5.refseq_genomic.gtf.gz" \
    --mane-ensembl-grch38-gtf "$workspace/gtf/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz" \
    --ensembl75-grch37-gtf "$workspace/gtf/Homo_sapiens.GRCh37.75.gtf.gz" \
    --grch38-refseq-fasta "$mane_refseq_genome" \
    --grch38-fasta "$mane_ensembl_genome" \
    --grch37-fasta "$ens75_genome" \
    --grch37-refseq-gtf "$workspace/gtf/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" \
    --refseq-output "$workspace/transcripts/MANE_1.5_GRCh37_REFSEQ_candidates.txt" \
    --ensembl-output "$workspace/transcripts/MANE_1.5_GRCh37_ENSEMBL_candidates.txt" \
    --report "$workspace/reports/MANE_1.5_GRCh37.mapping-review.tsv"

  drop_prepared "$workspace/fasta/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
  drop_prepared "$workspace/fasta/GRCh38.p14.genome.fa.gz"
  drop_prepared "$workspace/fasta/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz"

  gr37_refseq_out="$workspace/build/mane-1.5-grch37-refseq"
  gr37_refseq_run="$workspace/tmp/mane-1.5-grch37-refseq.run"
  mkdir -p "$gr37_refseq_out" "$gr37_refseq_run"
  cp -f "$workspace/gtf/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" "$gr37_refseq_out/"
  gr37_refseq_url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
  (
    cd "$gr37_refseq_run"
    PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}" \
      python3 "$repo_root/cava/RefSeqDB.py" --no_hg19 -b GRCh37 --nm-only \
      -i "$workspace/transcripts/MANE_1.5_GRCh37_REFSEQ_candidates.txt" \
      -D "$gr37_refseq_out" -o CAVA_MANE_1.5_GRCh37_REFSEQ.experimental \
      --url_gtf "$gr37_refseq_url"
  ) 2>&1 | tee "$workspace/logs/MANE_1.5_GRCh37_REFSEQ.experimental.log"
  gr37_refseq_base="$workspace/output/CAVA_MANE_1.5_GRCh37_REFSEQ.experimental"
  copy_products "$gr37_refseq_out/CAVA_MANE_1.5_GRCh37_REFSEQ.experimental.gz" \
    "$gr37_refseq_out/CAVA_MANE_1.5_GRCh37_REFSEQ.experimental.txt" "$gr37_refseq_base"
  create_seleno_refseq "$workspace/transcripts/MANE_1.5_GRCh37_REFSEQ_candidates.txt" \
    mane_1_5_grch37_refseq "$gr37_refseq_base.cesis"
  register_catalog mane-1.5-grch37-refseq-experimental "$gr37_refseq_base" experimental

  gr37_ens_build="$workspace/build/mane-1.5-grch37-ensembl"
  run_ensembl_builder 75 GRCh37 "$ens75_url" \
    "$workspace/gtf/Homo_sapiens.GRCh37.75.gtf.gz" \
    "$workspace/transcripts/MANE_1.5_GRCh37_ENSEMBL_candidates.txt" \
    CAVA_MANE_1.5_GRCh37_ENSEMBL.experimental "$gr37_ens_build"
  gr37_ens_base="$workspace/output/CAVA_MANE_1.5_GRCh37_ENSEMBL.experimental"
  copy_products "$gr37_ens_build/CAVA_MANE_1.5_GRCh37_ENSEMBL.experimental.gz" \
    "$gr37_ens_build/CAVA_MANE_1.5_GRCh37_ENSEMBL.experimental.txt" "$gr37_ens_base"
  create_seleno_ensembl "$ens75_map" mane_1_5_grch37_ensembl "$gr37_ens_base.cesis"
  register_catalog mane-1.5-grch37-ensembl-experimental "$gr37_ens_base" experimental
fi

python3 "$script_dir/validate_manifest.py" --repo "$repo_root" --check-files
cat <<EOF
Catalog build completed locally. Nothing was staged or pushed.

Review:
  $workspace/reports
  $workspace/logs

Publish selected catalog objects only with:
  scripts/catalogs/stage_and_push_catalogs.sh --no-push CATALOG_NAME [...]
EOF
