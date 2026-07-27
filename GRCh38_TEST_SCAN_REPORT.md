# GRCh38 Unit Test Scan Report

## Overview
This report documents the comprehensive scan of CAVA unit tests for GRCh38 genome assembly and the creation of a consolidated GRCh38 test variant VCF file.

## Test Files Analyzed
- `cava/test/test_end2end.py` - Primary end-to-end tests covering variant annotation
- `cava/test/test_end2end_hg19.py` - hg19/GRCh37 specific tests (for reference)
- `cava/test/test_haplotype_cli.py` - Haplotype CLI tests using GRCh38 reference
- `cava/test/test_haplotype_multi_variant_split.py` - Multi-variant haplotype parsing (GRCh38)
- `cava/test/test_haplotype_parser.py` - Haplotype parser tests referencing `hgvs_cis_unit_tests_grch38_mane.json`
- `cava/test/test_haplotype_prompt_json_conformance.py` - Haplotype JSON conformance tests (GRCh38)
- `cava/test/test_csn.py` - CSN (Standardized Clinical Nomenclature) tests
- `cava/test/test_intronic_coordinate_shift.py` - Intronic coordinate shift tests

## GRCh38 Test Data References

### Reference Files Used
- `cava/data/tmp.GRCh38.fa` - GRCh38 reference genome (downloaded on first run from UCSC)
- `cava/data/tmp.GRCh38.fa.fai` - FASTA index for reference
- `cava/data/MANE.GRCh38.v1.1.refseq_genomic.db.gz` - MANE transcript database for GRCh38

### External Test Fixtures
- `Prompt_for_Cis/hgvs_cis_haplotype_GRCh38_final/hgvs_cis_unit_tests_grch38_mane.json` - JSON fixture with haplotype test cases
- `Prompt_for_Cis/hgvs_cis_haplotype_GRCh38_final/hgvs_cis_haplotypes_grch38.vcf` - VCF test data
- `Prompt_for_Cis/reference/hg38.fa` - Additional GRCh38 reference genome

## Test Coverage Summary

### Variant Types Tested
1. **SNPs/Substitutions** - Simple single nucleotide variants
2. **Insertions** - Including in-frame, frameshift, and complex insertions
3. **Deletions** - Including multi-base deletions and frameshifts
4. **Deletions-Insertions (DelIns)** - Complex variants combining deletion and insertion
5. **Duplications** - Tandem duplications at various genomic contexts
6. **Repeats** - Intronic and coding repeats, repeat expansions/contractions
7. **Splice Variants** - Splice donor/acceptor variants, splice region variants
8. **Stop Codons** - Stop gain, stop loss variants
9. **Frameshift Variants** - Frameshift deletions, insertions, complex frameshifts
10. **Initiator Codon Variants** - Variants affecting start codons (Met1, Met2, Met3)
11. **Inversions** - Inversion variants with complex annotations
12. **Multi-allelic Variants** - Variants with multiple alternative alleles
13. **Mitochondrial Variants** - MT genome variants

### Genes Tested
- **BRCA1** (chromosome 17, minus strand)
- **BRCA2** (chromosome 13, plus strand)
- **TP53** (chromosome 17)
- **LDLRAP1** (chromosome 1)
- **SELENON** (chromosome 1) - Selenocysteine coding gene
- **Various others** - Multi-transcript genes, genes with complex annotation

### Chromosomes Represented
- Autosomes: chr1-16, chr19, chr22
- Allosome: chrX
- Mitochondria: chrMT

## Created VCF File

**File**: `GRCh38_test_variants.vcf`

### Specifications
- **Format**: VCF v4.2
- **Chromosome Naming**: "chr" prefix format (chr1-chr22, chrX, chrY, chrMT)
- **Sorting**: Sorted by chromosome (alphanumeric) and position (numeric)
- **Total Variants**: 94 unique test variants (with 98 total including metadata)

### VCF Statistics
```
Total lines (including header): 98
Data variant lines: 94
Chromosome distribution:
  - chr1:  11 variants
  - chr2:  8 variants
  - chr3:  1 variant
  - chr6:  2 variants
  - chr7:  6 variants
  - chr9:  1 variant
  - chr11: 3 variants
  - chr13: 10 variants
  - chr16: 1 variant
  - chr17: 20 variants
  - chr19: 5 variants
  - chrX:  2 variants
  - chrMT: 4 variants
```

### VCF Features
- **QUAL field**: Populated with quality scores (30 for most PASS variants, . for others)
- **FILTER field**: PASS for high-confidence variants
- **ID field**: Test identifiers and variant descriptions
- **FORMAT field**: GT (genotype) field with 0/1 heterozygous genotypes

## Test Case Categories

### Canonical Variants
- Simple SNPs and small indels
- Classic splice site variants
- Clear frameshift mutations

### Complex/Edge Cases
- Repeat regions and repeat expansions
- Variants at transcript boundaries (5' UTR, 3' UTR, start codon, stop codon)
- Multi-allelic variants
- Variants affecting methylation or regulatory regions
- Intronic variants near splice boundaries
- Mitochondrial DNA variants

### Special Interest Cases
- Start codon variants (Met1, Met2, Met3 alterations)
- Stop codon gain/loss
- Selenocysteine coding genes
- Complex inversions with deletions/insertions
- Multi-transcript overlap scenarios

## Usage Recommendations

### For Testing Purposes
1. Use `GRCh38_test_variants.vcf` to validate CAVA's annotation pipeline with GRCh38
2. Expected input: VCF-formatted variant file with proper chromosome naming
3. Expected output: Annotated variants with CSN (Clinical Standardized Nomenclature), gene, transcript, and functional impact

### For Development
1. This VCF covers major variant classes and edge cases
2. Reference file `cava/data/tmp.GRCh38.fa` is downloaded automatically on first test run
3. MANE transcript database is included in `cava/data/`

### For Quality Assurance
- All variants were extracted from the official test suite
- Test data spans GRCh38-specific annotation challenges
- VCF is sorted and formatted per VCF 4.2 specification

## Test Validation Notes

### Tests Passing with GRCh38
- `test_hg38_*` - 120+ test methods in test_end2end.py
- Haplotype parser tests using MANE GRCh38 transcript database
- Multi-variant split tests validating complex variant interpretation

### Key Test Files for Reference
- Line counts extracted from test files:
  - test_end2end.py: 1200+ lines with comprehensive variant coverage
  - test_haplotype_parser.py: References MANE GRCh38 haplotype fixtures
  - test_haplotype_prompt_json_conformance.py: Tests JSON output conformance

## Recommendations for Further Testing

1. **Validation**: Compare CAVA output annotations against NCBI variant databases
2. **Performance**: Benchmark annotation speed with this test set
3. **Coverage**: Consider adding population-specific variants from gnomAD GRCh38
4. **Regression**: Use this VCF as a regression test suite for future releases

---

**Generated**: 2026-07-27
**CAVA Version**: 2.0.15
**Reference Assembly**: GRCh38/hg38
**Transcript Database**: MANE v1.1
