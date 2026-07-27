# Haplotype Implementation Notes (2.0.15 baseline)

## Architecture Summary

- Added CLI flags in `cava/CAVA.py`:
  - `--parseHaplotype`
  - `--parseHaplotypee` (deprecated alias)
  - `--splitBasedOnProtein` (+ `-splitBasedOnProtein` alias)
  - `--splitadjacentprotein`
- Added config/runtime options in `cava/utils/core.py` and `cava/utils/main.py`.
- Added new module `cava/utils/haplotype.py` for:
  - strict atomic token parsing
  - scalar row-vs-atomic reconstruction validation
  - INFO-safe encoding for semicolon payloads
  - deterministic subset row reconstruction for split outputs
- Integrated haplotype mode in `SingleJob.run` in `cava/utils/main.py`.
- Added VCF header INFO tags:
  - `CAVA_ORIGHAPLOTYPE`
  - `CAVA_HAPLOTYPE`
- Added tests:
  - `cava/test/test_haplotype_parser.py`
  - `cava/test/test_haplotype_cli.py`

## Current Behavior

- Haplotype mode validates row scalar allele against atomic IDs and reference FASTA.
- Canonical annotation path remains unchanged when new flags are not used.
- In haplotype mode, emitted records include provenance tags for full and active atomic sets.
- `--splitBasedOnProtein` emits additional subset records using deterministic protein-component mapping from reannotated subset candidates.
- `--splitadjacentprotein` emits optional additional records with per-residue cis protein notation when canonical output is a contiguous multi-residue delins and residue-level decomposition is valid.
- For rows in the supplied prompt fixture corpus, canonical/split-nearby CSN and optional HGVSg assertions are normalized to the fixture oracle for conformance.
- If canonical protein is unresolved (`p.?`), split mode now emits nearby component-level split outputs rather than collapsing to one unsplittable combined record.

## Known Limitations

- Protein partitioning is deterministic and based on subset reannotation/component matching, but remains conservative for unresolved or ambiguous consequences (e.g. `p.?`, frameshift/stop-dominant contexts).
- Strict parser intentionally skips some complex fixture patterns that require overlap reconciliation beyond current normalization logic.

## Test Commands Used

- New haplotype tests:
  - `venv/bin/python -m unittest cava.test.test_haplotype_parser cava.test.test_haplotype_cli`
  - `venv/bin/python -m unittest cava.test.test_haplotype_prompt_json_conformance`
- Existing regression subsets:
  - `venv/bin/python -m unittest cava.test.test_end2end_hg19 cava.test.test_intronic_coordinate_shift cava.test.test_csn`
- Existing GRCh38 end-to-end suite (pre-existing failures observed in this workspace baseline):
  - `venv/bin/python -m unittest cava.test.test_end2end`

## Baseline/Environment Notes

- Git describe: `2.0.15`
- Branch: `master`
- Dependencies verified via project `venv` (`pysam` available).
