Post-Build Checklist

Purpose
- Ensure generated catalogs are consistent with companion transcript->protein maps.
- Ensure CAVA annotation remains stable when protein information is absent.
- Ensure CAVA annotation remains stable when CDS start/length are missing/invalid.

1. Artifact Integrity
- Verify for each catalog base name:
  - <name>.gz exists and is non-empty.
  - <name>.gz.tbi exists and is non-empty.
  - <name>.txt exists and is non-empty.
  - <name>.cesis exists and is non-empty.
- Run: gzip -t <name>.gz

Suggested command:
- for b in catalogs/files/*.gz; do gzip -t "$b" || exit 1; done

2. Catalog vs Transcript-Map Consistency
- Run validate_catalog_companions.py for each built catalog.
- PASS criteria:
  - missing_map_row = 0
  - missing_protein_id_with_cds = 0
- INFO-only (acceptable):
  - missing_protein_id_non_coding > 0 (non-coding transcripts with empty protein ID)

Suggested command:
- python3 scripts/catalogs/validate_catalog_companions.py --catalog catalogs/files/<name>.gz --transcript-map catalogs/files/<name>.txt --report catalog-build/reports/<name>.companion-map.tsv

3. CDS-Aware Missing Protein QC
- Confirm every transcript lacking protein ID in <name>.txt has no CDS in <name>.gz.
- If any transcript has CDS in catalog but no protein ID in map, mark as QC failure.

4. GSTT1 Exception Coverage
- Confirm GSTT1 appears in:
  - validated transcript lists
  - final catalog(s)
  - companion map
- Confirm contig/alias resolution supports the non-primary contig location without collapsing primary chromosome aliases.

Suggested commands:
- zgrep -n "\tGSTT1\t|^.*GSTT1.*$" catalogs/files/<name>.gz | head
- grep -n "GSTT1" catalogs/files/<name>.txt | head
- grep -n "GSTT1" catalog-build/validated-transcripts/*.txt | head

5. Runtime Safety: Missing Protein Mapping
- Run CAVA annotation on a variant in a transcript with no protein mapping entry.
- PASS criteria:
  - Process exits successfully.
  - HGVSp is rendered as '.' or equivalent safe placeholder.
  - Warning messages are acceptable; crashes/exceptions are not.

6. Runtime Safety: Missing/Invalid CDS Coordinates
- Build a small negative test using transcript entries with missing/invalid CDS positions.
- PASS criteria:
  - Annotation engine does not crash.
  - Transcript consequence logic degrades gracefully (e.g., non-coding or unknown protein impact markers).

Implementation note:
- Parser/runtime now treats missing CDS fields (e.g., '.', blank, NA) as non-coding and returns empty protein sequence instead of crashing.

7. Regression Spot Checks
- Compare transcript counts before and after filtering for each catalog.
- Review source and final sequence validation reports for unexpected spikes in failures.

8. External Coordinate Verification (Conditional)
- Use this check only in these situations:
  - during this build/debug phase,
  - when coordinate parsing/transformation code is changed,
  - or when a runtime discrepancy suggests coordinate mismatch.
- Verify a representative transcript sample against an external source (Ensembl REST).
- PASS criteria:
  - coordinate mismatches = 0
  - strand mismatches = 0
- Recommended command:
- python3 scripts/catalogs/verify_external_coordinates.py --catalog catalogs/files/CAVA_GENCODE_50_GRCh38.gz --report catalog-build/reports/gencode-50-grch38.external-coordinate-check.tsv --sample-size 50 --transcript ENST00000634102.1

9. Report and Sign-off
- Record per-catalog summary:
  - transcript rows in catalog
  - rows in map
  - missing_map_row
  - missing_protein_id_with_cds
  - missing_protein_id_non_coding
- Mark build as publishable only when all QC failures are zero.

Current policy summary:
- Empty protein ID is acceptable only when the transcript has no CDS in the catalog.
- Empty protein ID is a QC failure when the transcript has CDS in the catalog.
- GSTT1 must be retained and validated even when located on non-primary assembled sequence.
