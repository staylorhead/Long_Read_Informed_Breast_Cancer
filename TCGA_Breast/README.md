
# TCGA BRCA Breast Tumor Analysis Pipeline

This pipeline adapts the GTEx breast mammary tissue workflow for TCGA BRCA primary solid tumor samples. The overall structure is similar, but there are several TCGA-specific differences: genotypes come from imputed Affymetrix SNP6 arrays (rather than WGS), a tumor-specific long-read assembly (`veigaTum`) is used alongside GENCODE v45, and an extra BED reformatting step is included to handle TCGA's multi-sample-per-participant structure. Script numbering is shifted by ~1 relative to the GTEx pipelines to accommodate these differences.

---

## Preprocessing (`other_s0`–`other_s3`)

| Script | Description |
|--------|-------------|
| `other_s0_combine_gencodev45_tumorAssembly.R` | Combines GENCODE v45 with the Veiga ESPRESSO tumor long-read assembly GTF into a merged annotation. Assigns gene IDs to novel transcripts using SQANTI3 classification tables, with fallback to coordinate-based IDs for unannotated novels. |
| `other_s1_build_salmon_index_for_requantifying_TCGA_BRCA_fastqs.sh` | Builds a Salmon decoy-aware index for the combined GENCODE v45 + ESPRESSO tumor-filtered annotation using `gffread` + `generateDecoyTranscriptome.sh`. |
| `other_s2_run_salmon_quant.sh` | Salmon quantification of TCGA BRCA FASTQs against the combined index, using sequence bias correction (`--seqBias`). Sample IDs are resolved from a hash-code-to-TCGA-barcode lookup table. |
| `other_s2_run_salmon_quant_with_gcbias.sh` | Same as above with additional GC bias correction (`--gcBias`). Used for comparison or as the primary quantification depending on parameter settings. |
| `other_s3_make_se.R` | Loads Salmon output into a `SummarizedExperiment` via `tximeta`, linking to the combined tumor GTF. Also runs `catchSalmon` to capture overdispersion estimates. Saves the SE object for use in `s03`. |

**Sankey diagram scripts** (`other_s09`–`other_s12`): Build STAR indices and align TCGA tumor FASTQs against GENCODE v45 and the ESPRESSO tumor long-read assembly, extract read-to-transcript assignment tables, and compute read-mapping transition statistics for Sankey flow diagrams.

---

## Main Pipeline

### `s00` — Genotype QC

`s00_genotype_qc.sh` — Applies variant-level filters to TCGA BRCA imputed Affymetrix SNP6 BCF data with `plink2`: MAF > 1%, HWE p > 1e-6, genotype missingness < 5%. Extracts imputation quality scores (R², ER²) from the BCF INFO field. Produces two QC-filtered pgen sets: one with all passing variants and one additionally filtered to imputation R² > 0.8. Computes KING kinship coefficients and sample-level missingness for both.

### `s01` — Relatedness Check

`s01_relatedness_check.R` — Reads KING kinship output and identifies related sample pairs (kinship > 0.0884; 54 pairs found). Parses TCGA barcodes to extract participant IDs, sample types (tumor vs. normal), and TSS codes, and merges with subtype data for downstream use in sample selection.

### `s02` — Genotype VCF Preparation

`s02_make_genotype_vcfs.sh` — Computes genotype PCs on unrelated samples (LD pruning, then `plink2 --pca 50`), exports QC-filtered autosomal genotypes to bgzipped VCF with participant-level IDs (trimming full TCGA barcodes to 12-character participant IDs), tabix-indexes, and adds `chr` prefixes for QTLtools compatibility. Exports both GT-only and GT+DS (dosage) VCF versions.

### `s03` — Expression Preparation

`s03_prep_expr.R` — Loads the tximeta SE for a given annotation (`gencodev45`, `veigaTum`, or `comb`) and applies the normalization pipeline (ODS scaling, TMM, VST, length-scaled TPM as controlled by argument flags). Key TCGA-specific steps relative to the GTEx pipeline:
- Subsets to primary solid tumor samples only (sample type `01`)
- Applies ComBat batch correction on TSS (tissue source site), which is meaningful for TCGA but was disabled for GTEx
- Uses `TCGAbiolinks` for TCGA metadata integration
- Outputs separate `_tumor`-suffixed BED files

### `s04` — BED File Reformatting

`s04_reformat_bed.R` — TCGA-specific step with no GTEx equivalent. TCGA has multiple RNA-seq samples per participant; this script filters each expression BED to primary solid tumor samples only, then selects the single best-quality sample per participant based on alignment metrics (from `BRCA_txomeMetrics.txt`), producing a deduplicated, participant-level BED file for QTLtools.

### `s05` — Hidden Covariate Estimation (HCP)

`s05_HCP.sh` → **`child_HCP_by_param.R`** — Estimates hidden covariates (HCP or PCA-based) and combines them with genotype PCs to produce covariate files for QTLtools. Also bgzips/tabix-indexes expression BED files if not already done.

### `s06` — QTLtools Permutation Pass

`s06_QTLtools_perm_tx_by_chr.sh` / `s06_QTLtools_perm_gene_by_chr.sh` — Array job (22 chromosomes) running `QTLtools cis --permute` at transcript and gene level to map *cis*-eQTLs with permutation-based significance. No child R script; calls QTLtools directly.

### `s07` — Aggregate Permutation Results

`s07_agg_res_perm.sh` → **`child_qtltools_runFDR_cis.R`** — Concatenates per-chromosome permutation results and runs FDR correction (5% and 20%) to identify eGenes/eTx.

### `s08` — QTLtools Conditional Pass

`s08_QTLtools_conditional_tx_by_chr.sh` / `s08_QTLtools_conditional_gene_by_chr.sh` — Array job running `QTLtools cis --mapping` to identify independent secondary eQTL signals. No child R script.

### `s09` — Aggregate Conditional Results

`s09_agg_res_conditional.sh` — Concatenates per-chromosome conditional results into a single file per level. No child R script.

### `s10` — Fine-mapping (Univariate SuSiE)

`s10_finemapping_univariate.sh` → **`child_finemapping_univariate.R`** — Array job running univariate SuSiE fine-mapping on eQTL loci to generate PIPs and credible sets.

### `s11` — Colocalization with GWAS

Runs nominal QTLtools per chunk (1 of 100) and pipes summary statistics into colocalization with BC GWAS. Nominal eQTL results are deleted after running colocalization to save space.

| Script | Child script | GWAS phenotypes | Level |
|--------|-------------|----------------|-------|
| `s11_coloc_BC_gene.sh` | `child_run_coloc_BC.R` | overallBC, lumA, lumB, lumBHER2neg, HER2enriched, TN | gene only |
| `s11_coloc_BC_tx.sh` | `child_run_coloc_BC.R` | overallBC, lumA, lumB, lumBHER2neg, HER2enriched, TN | tx only |

### `s12` — isoTWAS Model Training

`s12_isoTWAS_model_overlap_ss_munged.sh` → **`child_fit_isoTWAS_model_overlap_ss_munged.R`** — Array job fitting isoTWAS predictive models of transcript expression from genotype, called per gene bin.

### `s13` — isoTWAS Association Testing

`s13_eval_isoTWAS_BC_overlap_ss_munged.sh` → **`child_eval_isoTWAS_BC_overlap_ss_munged.R`** — Array job running isoTWAS association testing against BC/MD GWAS using precomputed LD and models from `s12`.

### `s14` — isoTWAS Screen-and-Confirm

`s14_isoTWAS_ScreenConfirm_overlap_ss_munged.R` — Self-contained R script. Applies the screen-and-confirm multiple testing procedure across all traits and quantification methods (`gencodev45`, `veigaTum`, `comb`). Filters the HLA region (chr6:27–35Mb), computes gene-level screen p-values, FDR-adjusts, and derives isoform-level confirmation p-values. Writes full results and significant associations to TSV files.

### `s15` — Aggregate isoTWAS Significant Hits (no fine-mapping filter)

`s15_isoTWAS_noFM_overlap_ss_munged.R` — Self-contained R script. Collects `SignificantAssociations_isoTWAS_noFineMap.tsv` files across all quantification methods and traits into a single combined output table.

### `s16` — eQTL LD Computation (for TWAS input)

`s16_eqtl_ld_gwas.sh` → **`child_eqtl_ld_gwas.R`** — Array job computing LD matrices between eQTL variants and GWAS SNPs for use in isoTWAS.

### `s17` — Aggregate Fine-mapping Results

`s17_agg_res_finemap.R` — Self-contained R script. Reads per-locus SuSiE `.RData` files, extracts credible sets and high-PIP SNPs, and saves aggregated gene- and tx-level summaries.

### `s18` — Fine-mapping LD Computation (for GWAS overlap)

`s18_finemap_ld_gwas.sh` → **`child_finemap_ld_gwas.R`** — Computes LD between fine-mapped eQTL credible set SNPs and GWAS variants for colocalization support.

### `s19` — Aggregate Colocalization Results

`s19_agg_res_coloc.R` — Self-contained R script. Reads all per-chunk coloc `.RData` files, extracts PP.H0–H4 and top colocalizing SNPs (by SNP.PP.H4), and combines into a single aggregated results table.

### `s20` — Standalone Nominal eQTL Pass (stored output)

`s20_QTLtools_nominal.sh` — Runs a nominal QTLtools pass at gene and tx level with output retained for downstream use. No child R script.
