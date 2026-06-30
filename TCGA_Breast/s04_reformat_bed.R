#!/usr/bin/env Rscript

########################################################################
# Set up R library path
########################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3", myPaths)
.libPaths(myPaths)

########################################################################
# Load dependencies
########################################################################
library(data.table)
library(dplyr)
library(SummarizedExperiment)

########################################################################
# Load tximeta object and extract sample metadata
########################################################################
load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/se_gencodev38.RData")
sample_ids_se <- colnames(assay(se))
sample_dat <- data.frame(barcode = sample_ids_se)

sample_dat$TSS <- sapply(strsplit(sample_dat$barcode, "-"), `[[`, 2)
sample_dat$participant <- sapply(strsplit(sample_dat$barcode, "-"), `[[`, 3)
sample_dat$sample_type <- substr(sapply(strsplit(sample_dat$barcode, "-"), `[[`, 4), 1, 2)

sample_dat$sample_type <- factor(sample_dat$sample_type,
  levels = c("01", "06", "11"),
  labels = c("Primary Solid Tumor", "Metastatic", "Solid Tissue Normal"))
sample_dat$barcode_short <- substr(sample_dat$barcode, 1, 12)

# Clean up
rm(se)

########################################################################
# Get tumor sample IDs
########################################################################
tumor_ids <- sample_dat$barcode[sample_dat$sample_type == "Primary Solid Tumor"]

########################################################################
# Load gene expression files and metrics
########################################################################

file_list <- list.files("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed")
file_list <- file_list[grep("TCGA",file_list)]
file_list <- file_list[grep(".bed.gz",file_list)]

gene_file_list <- file_list[grep("^TCGA_gene", file_list)]
tx_file_list <- file_list[grep("^TCGA_tx", file_list)]

metrics <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/BRCA_txomeMetrics.txt")
metrics$prior_ids_sub <- sapply(strsplit(metrics$sample_ID, ".rna_seq.transcriptome.gdc_realn.bam"), `[[`, 1)

barcode_dat <- fread('/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/fastq_hash_TCGA_ID_Full_ID.txt')
names(barcode_dat) <- c("hash", "full_barcode", "short_barcode")

metrics <- merge(x = metrics, y = barcode_dat, by.x = "prior_ids_sub", by.y = "hash")

########################################################################
# Loop over gene bed files and filter to best tumor sample per participant
########################################################################

for (i in seq_along(gene_file_list)) {

  print(gene_file_list[i])
  
  # Load expression matrix
  bed_gene <- fread(paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/", gene_file_list[i]))
  print(ncol(bed_gene)-6)

  # Filter to tumor samples present in this matrix
  sample_match <- colnames(bed_gene) %in% tumor_ids
  bed_gene_tumor <- bed_gene[, c(1:6, which(sample_match)), with = FALSE]
  print(ncol(bed_gene_tumor)-6)

  # Extract and match sample names to participant IDs
  expr_cols <- names(bed_gene_tumor)[-1:-6]
  sample_dt <- data.table(sample_name = expr_cols)
  sample_dt[, participant_id := substr(sample_name, 1, 12)]

  metrics_dt <- as.data.table(metrics)
  sample_dt <- merge(sample_dt, metrics_dt, by.x = "sample_name", by.y = "full_barcode", all.x = TRUE)

  # Choose best sample per participant using highest PCT_USABLE_BASES
  repeat_ids <- names(table(sample_dt$participant_id))[table(sample_dt$participant_id)>1]
  non_repeat_ids <- names(table(sample_dt$participant_id))[table(sample_dt$participant_id)==1]
  best_samples <- sample_dt[sample_dt$participant_id %in% non_repeat_ids,]

  for(ii in 1:length(repeat_ids)){
    tmp <- sample_dt[sample_dt$participant_id==repeat_ids[ii]]
    keep <- which.max(tmp$PCT_USABLE_BASES)
    best_samples <- rbind(best_samples,tmp[keep,])
  }

 # setorder(sample_dt, participant_id, -PCT_USABLE_BASES)
 # best_samples <- sample_dt[, .SD[1], by = participant_id]

  keep_sample_names <- best_samples$sample_name
  fixed_cols <- names(bed_gene_tumor)[1:6]
  cols_to_keep <- c(fixed_cols, keep_sample_names)

  bed_gene_tumor_best <- bed_gene_tumor[, ..cols_to_keep]

  # Create output prefix based on file name (was missing in original script)
  prefix <- sub(".bed.gz", "", gene_file_list[i])

  print(ncol(bed_gene_tumor_best)-6)

    # Save full barcodes before renaming
  writeLines(
    keep_sample_names,
    con = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/", prefix, "_tumor_full_barcodes.txt")
  )

  setnames(
  bed_gene_tumor_best,
  old = keep_sample_names,
  new = substr(keep_sample_names, 1, 12)
  )

  # Save filtered matrix
  fwrite(bed_gene_tumor_best,paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/',prefix,'_tumor.bed'),
       sep='\t',col.names=T,row.names=F)
}

########################################################################
# Loop over tx bed files and filter to best tumor sample per participant
########################################################################

for (i in seq_along(tx_file_list)) {

  print(tx_file_list[i])
  
  # Load expression matrix
  bed_gene <- fread(paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/", tx_file_list[i]))
  print(ncol(bed_gene)-6)

  # Filter to tumor samples present in this matrix
  sample_match <- colnames(bed_gene) %in% tumor_ids
  bed_gene_tumor <- bed_gene[, c(1:6, which(sample_match)), with = FALSE]
  print(ncol(bed_gene_tumor)-6)

  # Extract and match sample names to participant IDs
  expr_cols <- names(bed_gene_tumor)[-1:-6]
  sample_dt <- data.table(sample_name = expr_cols)
  sample_dt[, participant_id := substr(sample_name, 1, 12)]

  metrics_dt <- as.data.table(metrics)
  sample_dt <- merge(sample_dt, metrics_dt, by.x = "sample_name", by.y = "full_barcode", all.x = TRUE)

  # Choose best sample per participant using highest PCT_USABLE_BASES
  setorder(sample_dt, participant_id, -PCT_USABLE_BASES)
  best_samples <- sample_dt[, .SD[1], by = participant_id]

  keep_sample_names <- best_samples$sample_name
  fixed_cols <- names(bed_gene_tumor)[1:6]
  cols_to_keep <- c(fixed_cols, keep_sample_names)

  bed_gene_tumor_best <- bed_gene_tumor[, ..cols_to_keep]

  # Create output prefix based on file name (was missing in your original script)
  prefix <- sub(".bed.gz", "", tx_file_list[i])

  print(ncol(bed_gene_tumor_best)-6)

      # Save full barcodes before renaming
  writeLines(
    keep_sample_names,
    con = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/", prefix, "_tumor_full_barcodes.txt")
  )

  setnames(
  bed_gene_tumor_best,
  old = keep_sample_names,
  new = substr(keep_sample_names, 1, 12)
  )

  # Save filtered matrix
  fwrite(bed_gene_tumor_best,paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/',prefix,'_tumor.bed'),
       sep='\t',col.names=T,row.names=F)
}