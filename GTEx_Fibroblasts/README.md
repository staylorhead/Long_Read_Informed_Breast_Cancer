# GTEx v8 Cultured Fibroblasts Analysis Pipeline

This pipeline analyzes GTEx v8 cultured fibroblast samples. The differences between this pipeline and that for healthy breast mammary tissue are the tissue type and a fibroblast-specific long-read assembly from GTEx v9 SQANTI3 (referred to as `ESPRESSO` rather than `veigaNorm`). Scripts are numbered in execution order; `other_` prefixed scripts are supporting/preprocessing steps.

---

## Preprocessing (`other_s00`–`other_s08`)

| Script | Description |
|--------|-------------|
| `other_s00_generate_bam_file_lists.R` | Filters GTEx sample attribute metadata to cultured fibroblast RNASEQ samples, writes a sample manifest used by downstream steps. |
| `other_s01_convert_bam_fastq.sh` | BAM → FASTQ conversion (array job, 504 samples). Sorts BAMs by read name with `samtools sort -n`, then extracts paired-end FASTQs with `samtools fastq`. |
| `other_s02_combine_gencodev45_novelAssembly.R` | Combines GENCODE v45 with the GTEx v9 SQANTI3 ESPRESSO fibroblast long-read assembly GTF into a single merged annotation. Handles novel transcript gene ID assignment using associated gene lookups from the SQANTI3 classification table. |
| `other_s03_build_salmon_index_for_requantifying_GTEx_fastqs_lr.sh` | Builds a Salmon decoy-aware index for the ESPRESSO fibroblast-filtered annotation using `gffread` + `generateDecoyTranscriptome.sh`. |
| `other_s03_build_salmon_index_for_requantifying_GTEx_fastqs_comb.sh` | Builds a Salmon decoy-aware index for the combined GENCODE v45 + ESPRESSO fibroblast annotation. |
| `other_s04_quant_fastq_gencodev45.sh` | Salmon quantification of FASTQs against the GENCODE v45 index (array job). |
| `other_s04_quant_fastq_lr.sh` | Salmon quantification against the ESPRESSO fibroblast-filtered index (array job). |
| `other_s04_quant_fastq_combined.sh` | Salmon quantification against the combined GENCODE v45 + ESPRESSO fibroblast index. |
| `other_s05_make_se.R` | Loads Salmon output into `SummarizedExperiment` objects via `tximeta`, linking to the appropriate GTF metadata. Saves `.RData` files for use in `s03`. |
| `other_s06_align.sh` | STAR alignment of trimmed FASTQs to GRCh38 (array job). Used for read-mapping diagnostics. |
| `other_s07_picard_tools.sh` | Runs Picard `CollectGcBiasMetrics` and related QC tools on STAR-aligned BAMs. |
| `other_s08_multiqc.sh` | Aggregates Picard QC output into a MultiQC report. |

**Sankey diagram scripts** (`other_s09`–`other_s12`): Build STAR indices and align fibroblast FASTQs against GENCODE v45, the ESPRESSO long-read assembly, and the combined annotation, then extract read-to-transcript assignment tables and compute read-mapping transition statistics. Includes a separate combined-annotation Sankey track (`other_s09_sankey_build_star_index_comb_fibro.sh`, `other_s10_sankey_align_comb_fibro.sh`, `other_s12_read_mapping_transitions_comb_fibro.R`).

---

## Main Pipeline

### `s03` — Expression Preparation

`s03_prep_expr.R` — Identical logic to the breast pipeline. Loads the tximeta `SummarizedExperiment` for a given annotation (`gencodev45`, `ESPRESSO`, or `comb`), applies the normalization pipeline (ODS scaling, TMM, VST, length-scaled TPM as controlled by argument flags), runs WGCNA-based outlier detection, and writes transcript- and gene-level BED files for QTLtools. Note: the long-read annotation label here is `ESPRESSO` (vs. `veigaNorm` in breast).

### `s04` — Hidden Covariate Estimation (HCP)

`s04_HCP.sh` → **`child_HCP_by_param.R`** — Estimates hidden covariates (HCP or PCA-based) and combines them with genotype PCs to produce covariate files for QTLtools. Also bgzips/tabix-indexes expression BED files if not already done.

### `s05` — QTLtools Permutation Pass

`s05_QTLtools_perm_tx_by_chr.sh` / `s05_QTLtools_perm_gene_by_chr.sh` — Array job (22 chromosomes) running `QTLtools cis --permute` at transcript and gene level to map *cis*-eQTLs with permutation-based significance. No child R script; calls QTLtools directly.

### `s06` — Aggregate Permutation Results

`s06_agg_res_perm.sh` → **`child_qtltools_runFDR_cis.R`** — Concatenates per-chromosome permutation results, then calls the FDR script twice per level (5% and 20% thresholds) to identify eGenes/eTx.

### `s07` — QTLtools Conditional Pass

`s07_QTLtools_conditional_tx_by_chr.sh` / `s07_QTLtools_conditional_gene_by_chr.sh` — Array job running `QTLtools cis --mapping` to identify independent secondary eQTL signals. No child R script.

### `s08` — Aggregate Conditional Results

`s08_agg_res_conditional.sh` — Concatenates per-chromosome conditional results into a single file per level. No child R script.

### `s09` — Fine-mapping (Univariate SuSiE)

`s09_finemapping_univariate.sh` → **`child_finemapping_univariate.R`** — Array job running univariate SuSiE fine-mapping on eQTL loci to generate PIPs and credible sets.

### `s10` — Colocalization with GWAS

Runs nominal QTLtools per chunk (1 of 100) and pipes summary statistics into colocalization with BC GWAS. Nominal eQTL results are deleted after running colocalization.

| Script | Child script | GWAS phenotypes | Level |
|--------|-------------|----------------|-------|
| `s10_coloc_BC_gene.sh` | `child_run_coloc_BC.R` | overallBC, lumA, lumB, lumBHER2neg, HER2enriched, TN | gene only |
| `s10_coloc_BC_tx.sh` | `child_run_coloc_BC.R` | overallBC, lumA, lumB, lumBHER2neg, HER2enriched, TN | tx only |

### `s11` — isoTWAS Model Training

`s11_isoTWAS_model_overlap_ss_munged.sh` → **`child_fit_isoTWAS_model_overlap_ss_munged.R`** — Array job fitting isoTWAS predictive models of transcript expression from genotype, called per gene bin.

### `s12` — isoTWAS Association Testing

`s12_eval_isoTWAS_BC_overlap_ss_munged.sh` → **`child_eval_isoTWAS_BC_overlap_ss_munged.R`** — Array job running isoTWAS association testing against BC/MD GWAS, using precomputed LD and the models from `s11`.

### `s13` — isoTWAS Screen-and-Confirm

`s13_isoTWAS_ScreenConfirm_overlap_ss_munged.R` — Self-contained R script. Applies the screen-and-confirm multiple testing procedure across all traits and quantification methods (`gencodev45`, `ESPRESSO`, `comb`). Filters the HLA region (chr6:27–35Mb), computes gene-level screen p-values, FDR-adjusts, and derives isoform-level confirmation p-values. Writes full results and significant associations to TSV files.

### `s14` — Aggregate isoTWAS Significant Hits (no fine-mapping filter)

`s14_isoTWAS_noFM_overlap_ss_munged.R` — Self-contained R script. Collects `SignificantAssociations_isoTWAS_noFineMap.tsv` files across all quantification methods and traits into a single combined output table.

### `s15` — eQTL LD Computation (for TWAS input)

`s15_eqtl_ld_gwas.sh` → **`child_eqtl_ld_gwas.R`** — Array job computing LD matrices between eQTL variants and GWAS SNPs for use in isoTWAS.

### `s16` — Aggregate Fine-mapping Results

`s16_agg_res_finemap.R` — Self-contained R script. Reads per-locus SuSiE `.RData` files for selected parameter combinations, extracts credible sets and high-PIP SNPs, and saves aggregated summaries for gene- and tx-level results separately.

### `s17` — Fine-mapping LD Computation (for GWAS overlap)

`s17_finemap_ld_gwas.sh` → **`child_finemap_ld_gwas.R`** — Computes LD between fine-mapped eQTL credible set SNPs and GWAS variants for colocalization support.

### `s18` — Aggregate Colocalization Results

`s18_agg_res_coloc.R` — Self-contained R script. Reads all per-chunk coloc `.RData` files, extracts PP.H0–H4 summary statistics and top colocalizing SNPs (by SNP.PP.H4), and combines into a single aggregated results table.

### `s19` — Standalone Nominal eQTL Pass (stored output)

`s19_QTLtools_nominal.sh` — Runs a nominal QTLtools pass at gene and tx level with output retained. Used for downstream analyses requiring full nominal summary statistics. No child R script.

