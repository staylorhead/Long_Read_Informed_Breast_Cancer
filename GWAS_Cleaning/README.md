
# GWAS Summary Statistics Munging Pipeline 

This pipeline formats raw BCAC 2020 breast cancer GWAS summary statistics with `MungeSumstats` and lifts them from GRCh37 to GRCh38. It runs as two stages on the HPC (LSF/`bsub`).

## Stage 1: Format raw summary statistics (`s00`)

**Scripts:** `s00_format_gwas_BC.sh` → `child_format_gwas_ss_BC.R`

- `s00_format_gwas_BC.sh` is the LSF wrapper: loads the R module and calls the formatting R script with `${build} ${pheno}`.
- `child_format_gwas_ss_BC.R` does the actual work, branching on phenotype:
  - **`overallBC`**: reads the iCOGS/OncoArray meta-analysis overall breast cancer summary stats, selects relevant columns, takes the minimum imputation R² across iCOGS/OncoArray, renames columns to `MungeSumstats`-expected format, and runs `format_sumstats()` with `INFO_filter = 0.3`, strand-ambiguous SNP removal, and biallelic filtering turned **off**.
  - **`lumA` / `lumB` / `lumBHER2neg` / `HER2enriched`**: reads the intrinsic-subtype meta-analysis file, pulls the subtype-specific beta/SE columns, computes a Wald-test p-value from beta/SE, and runs the same `format_sumstats()` call.
  - **`TN`** (triple negative): reads the CIMBA BRCA1/BCAC triple-negative meta-analysis file, renames columns, and runs the same `format_sumstats()` call.
- Output: `munged_<build>_<pheno>_ambig_strand_removed_no_biallelic_filt.tsv.gz`, written to the BCAC2020 data directory.

## Stage 2: Liftover to GRCh38 (`s01`)

**Scripts:** `s01_liftover_gwas_BC.sh` → `child_liftover_gwas_ss_BC.R`

- `s01_liftover_gwas_BC.sh` is the LSF wrapper: loads R and runs the liftover R script with `${ref} ${target} ${pheno}`.
- `child_liftover_gwas_ss_BC.R` reads the munged GRCh37 file(s) for a phenotype, runs `MungeSumstats::liftover()` (Ensembl chain files, NCBI chromosome style) to convert coordinates to GRCh38, and writes lifted output for both the full and strand-ambiguous-removed versions.

## Notes

- The liftover script as written expects munged files *without* the `_no_biallelic_filt` suffix (e.g. `munged_GRCh37_<pheno>.tsv.gz` and `..._ambig_strand_removed.tsv.gz`), while the current `s00` formatting step only writes the `_ambig_strand_removed_no_biallelic_filt` version — filenames will need to be reconciled (or the liftover script updated) for `s01` to find its input.


