# GTEx v8 Breast Mammary Tissue Analysis Pipeline

This pipeline processes GTEx v8 breast mammary tissue data from raw BAMs through eQTL mapping, fine-mapping, colocalization with breast cancer GWAS, and isoTWAS analysis. Scripts are numbered in execution order; `other_` prefixed scripts are supporting/preprocessing steps that don't fit neatly into the main sequence.

---

## Preprocessing (`other_s00`–`other_s08`)

These steps go from raw GTEx BAMs to Salmon quantification files and aligned BAMs, and are largely run once.

| Script | Description |
|--------|-------------|
| `other_s00_generate_bam_file_lists.R` | Filters GTEx sample attribute metadata to breast mammary tissue RNASEQ samples, writes a sample manifest used by downstream steps. |
| `other_s01_convert_bam_fastq.sh` | BAM → FASTQ conversion (array job, 440 samples). Sorts BAMs by read name with `samtools sort -n`, then extracts paired-end FASTQs with `samtools fastq`. |
| `other_s01_convert_bam_fastq_again_for_sankey.sh` | Repeat BAM → FASTQ conversion for a read-mapping Sankey diagram. |
| `other_s02_combine_gencodev45_normalAssembly.R` | Combines GENCODE v45 and the Veiga ESPRESSO long-read normal assembly GTFs into a single merged annotation for use as a custom transcriptome reference. |
| `other_s03_build_salmon_index_for_requantifying_GTEx_fastqs.sh` | Builds Salmon decoy-aware indices for the combined GTF annotation using `gffread` + `generateDecoyTranscriptome.sh`. |
| `other_s04_quant_fastq_gencodev45.sh` | Salmon quantification of FASTQs against the GENCODE v45 index (array job). |
| `other_s04_quant_fastq_veigaNorm.sh` | Salmon quantification against the Veiga ESPRESSO normal-filtered index (array job). |
| `other_s04_quant_fastq_combined.sh` | Salmon quantification against the combined GENCODE v45 + Veiga index. |
| `other_s05_make_se.R` | Loads Salmon output into `SummarizedExperiment` objects via `tximeta`, linking to GENCODE v45 or Veiga GTF metadata. Saves `.RData` files for use in `s03`. |
| `other_s06_align.sh` | STAR alignment of trimmed FASTQs to GRCh38 (array job). Runs Trim Galore first, then STAR with GENCODE v45 annotation. Used for read-mapping diagnostics. |
| `other_s07_picard_tools.sh` | Runs Picard `CollectGcBiasMetrics` and related QC tools on STAR-aligned BAMs. |
| `other_s08_multiqc.sh` | Aggregates Picard QC output into a MultiQC report. |

**Sankey diagram scripts** (`other_s09`–`other_s12`): Build STAR indices and align breast tissue FASTQs (both GTEx/GENCODE and long-read/Veiga), extract read-to-transcript assignment tables, and compute read-mapping transition statistics for Sankey flow diagrams comparing annotations. The R script `other_s12_read_mapping_transitions_breast.R` does the core statistical summarization.

---

## Main Pipeline

### `s00` — Genotype QC

`s00_genotype_qc.sh` — Applies variant-level filters to GTEx WGS data with `plink2`: MAF > 1%, genotype missingness < 5%, HWE p > 1e-6. Computes KING kinship coefficients and sample-level missingness.

### `s01` — Relatedness Check

`s01_relatedness_check.R` — Reads KING kinship output and flags related samples (kinship > 0.0884). No related samples were found in the GTEx breast cohort.

### `s02` — Genotype VCF Preparation

`s02_make_genotype_vcfs.sh` — Exports QC-filtered genotypes to bgzipped VCF (autosomes only), sorts and tabix-indexes, and adds `chr` prefixes to chromosome labels for compatibility with QTLtools.

### `s03` — Expression Preparation

`s03_prep_expr.R` — Loads the tximeta `SummarizedExperiment` for a given annotation (`gencodev45`, `veigaNorm`, or `comb`), applies a normalization pipeline controlled by flags passed as arguments:

- Remove EXCLUDE-flagged samples and chrX/Y/M transcripts
- Optionally compute length-scaled TPM (`lstpm`)
- Filter lowly expressed transcripts/genes (TPM > 0.1 in > 25% of samples)
- Optional ODS scaling (divide counts by overdispersion), TMM normalization, VST transformation
- WGCNA-based outlier detection and removal
- PCA of processed expression
- Writes transcript- and gene-level BED files for QTLtools input

### `s04` — Hidden Covariate Estimation (HCP)

`s04_HCP.sh` → **`child_HCP_by_param.R`** — Calls the child script to estimate hidden covariates (HCP or PCA-based, controlled by the `method` and `k` parameters), combining them with genotype PCs to produce covariate files for QTLtools. Also bgzips/tabix-indexes the expression BED files if not already done.

### `s05` — QTLtools Permutation Pass

`s05_QTLtools_perm_tx_by_chr.sh` / `s05_QTLtools_perm_gene_by_chr.sh` — Array job (22 chromosomes) running `QTLtools cis --permute` at transcript and gene level, respectively, to map *cis*-eQTLs with permutation-based significance. No child R script; calls QTLtools directly.

### `s06` — Aggregate Permutation Results

`s06_agg_res_perm.sh` → **`child_qtltools_runFDR_cis.R`** — Concatenates per-chromosome permutation results, then calls the FDR script twice per level (5% and 20% thresholds) to identify eGenes/eTx.

### `s07` — QTLtools Conditional Pass

`s07_QTLtools_conditional_tx_by_chr.sh` / `s07_QTLtools_conditional_gene_by_chr.sh` — Array job running `QTLtools cis --mapping` to identify independent secondary eQTL signals conditional on the top eQTL per feature. No child R script.

### `s08` — Aggregate Conditional Results

`s08_agg_res_conditional.sh` — Concatenates per-chromosome conditional results into a single file for each level (gene, tx). No child R script.

### `s09` — Fine-mapping (Univariate SuSiE)

`s09_finemapping_univariate.sh` → **`child_finemapping_univariate.R`** — Array job calling the child R script, which runs univariate SuSiE fine-mapping on eQTL loci, generating posterior inclusion probabilities (PIPs) and credible sets.

### `s10` — Colocalization with GWAS

`s10_coloc_BC.sh` → **`child_run_coloc_BC.R`** — Run a nominal QTLtools pass per chunk and immediately pipe the summary statistics into colocalization with GWAS. Nominal eQTL results are deleted after coloc to save disk space.

### `s11` — isoTWAS Model Training

`s11_isoTWAS_model_overlap_ss_munged.sh` → **`child_fit_isoTWAS_model_overlap_ss_munged.R`** — Array job fitting isoTWAS predictive models of transcript expression from genotype, using the processed expression and covariate files. The child script is called per gene bin.

### `s12` — isoTWAS Association Testing

`s12_eval_isoTWAS_BC_overlap_ss_munged.sh` → **`child_eval_isoTWAS_BC_overlap_ss_munged.R`** — Array job running isoTWAS association testing for a given GWAS (BC or MD), using precomputed LD and the isoTWAS models from `s11`.

### `s13` — isoTWAS Screen-and-Confirm

`s13_isoTWAS_ScreenConfirm_overlap_ss_munged.R` — Self-contained R script (no child). Applies the isoTWAS screen-and-confirm multiple testing procedure across all traits and quantification methods (`gencodev45`, `veigaNorm`, `comb`). Filters the HLA region (chr6:27–35Mb), computes gene-level screen p-values, FDR-adjusts, and derives isoform-level confirmation p-values. Writes full results and significant associations to TSV files.

### `s14` — Aggregate isoTWAS Significant Hits (no fine-mapping filter)

`s14_isoTWAS_noFM_overlap_ss_munged.R` — Self-contained R script. Collects the `SignificantAssociations_isoTWAS_noFineMap.tsv` files from all quantification methods and traits into a single combined output table.

### `s15` — eQTL LD Computation (for TWAS input)

`s15_eqtl_ld_gwas.sh` → **`child_eqtl_ld_gwas.R`** — Array job computing LD matrices between eQTL variants and GWAS SNPs for use in isoTWAS. The child script is called per chromosome.

### `s16` — Aggregate Fine-mapping Results

`s16_agg_res_finemap.R` — Self-contained R script. Reads per-locus SuSiE fine-mapping `.RData` files for selected parameter combinations, extracts credible sets and high-PIP SNPs, and saves aggregated summaries for gene- and tx-level results separately.

### `s17` — Fine-mapping LD Computation (for GWAS overlap)

`s17_finemap_ld_gwas.sh` → **`child_finemap_ld_gwas.R`** — Similar to `s15`, computes LD between fine-mapped eQTL credible set SNPs and GWAS variants for colocalization support.

### `s18` — Aggregate Colocalization Results

`s18_agg_res_coloc.R` — Self-contained R script. Reads all per-chunk coloc `.RData` files, extracts summary statistics (PP.H0–H4) and top colocalizing SNPs (by SNP.PP.H4), and combines into a single aggregated results table.

### `s19` — Standalone Nominal eQTL Pass (stored output)

`s19_QTLtools_nominal.sh` — Runs a nominal QTLtools pass at gene and tx level and retains the output (unlike `s10`, which deletes nominal results after coloc). No child R script; calls QTLtools directly. Used for downstream analyses requiring full nominal summary statistics.
