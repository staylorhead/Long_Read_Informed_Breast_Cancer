#!/usr/bin/env Rscript

########################################################################
# change library to local and load dependencies
########################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

library(data.table)
library(LDlinkR)
library(igraph)
library(rtracklayer)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
assembly <- as.character(args[1])
level <- as.character(args[2])
k <- as.numeric(args[3])
finemap_dir <- as.character(args[4]) 
chr <- as.numeric(args[5])
L <- as.numeric(args[6])

# assembly="gencodev45"
# level="gene"
# k=15
# finemap_dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap"
# chr=1
# L=10

####################################################################################
# begin code
####################################################################################

gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_tumor_filtered.gtf"
gtf <- import(gtf_file)
df <- data.frame(
  seqnames = as.character(seqnames(gtf)),
  transcript_id = mcols(gtf)$transcript_id,
  gene_id = mcols(gtf)$gene_id
)
df <- unique(df)
if (level == "gene") {
  # Subset to unique genes and their chromosomes
  df_unique <- unique(df[, c("gene_id", "seqnames")])
} else if (level == "tx") {
  # Keep all transcript-level rows as-is
  df_unique <- df
  df_unique <- df_unique[,c("transcript_id","seqnames")]
  colnames(df_unique)[1] <- "gene_id"
} else {
  stop("Unexpected level value. Expected 'gene' or 'tx'.")
}

# read in finemapping results
load(paste0(finemap_dir,"/TCGA_aggregated_finemap_res_byL_",level,".RData"))
summary_results$level <- unlist(lapply(strsplit(summary_results$set_name,"_",), `[[`, 2))
summary_results$annot <- unlist(lapply(strsplit(summary_results$set_name,"_",), `[[`, 3))
summary_results$L <- unlist(lapply(strsplit(summary_results$set_name,"_",), `[[`, 10))
summary_results <- summary_results[summary_results$annot==assembly,]
filt <- paste0("L",L)
summary_results <- summary_results[summary_results$L==filt,]

# Now merge with summary_results by gene_id or transcript_id accordingly
if (level == "gene") {
  merged_df <- merge(summary_results, df_unique, by.x = "transcript_id", by.y = "gene_id", all.x = TRUE)
} else if (level == "tx") {
  merged_df <- merge(summary_results, df_unique, by.x = "transcript_id", by.y = "gene_id", all.x = TRUE)
}
dat <- merged_df[merged_df$seqnames==paste0("chr",chr),]

# Load and subset GWAS loci
dat_bc <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/gwas_loci/indep_BC_gwas_loci_332.txt"))
dat_bc <- dat_bc[which(dat_bc$CHR == chr), ]

if (nrow(dat_bc) == 0) {
  stop("No GWAS loci on this chromosome — exiting job.")
}

pvar <- data.frame(fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr,".bim")))
# pvar$chr_pos <- paste0("chr",pvar$V1,":",pvar$V2)
pvar$id1 <- paste(pvar$V1,pvar$V4,pvar$V5,pvar$V6,sep=":")
pvar$id2 <- paste(pvar$V1,pvar$V4,pvar$V6,pvar$V5,sep=":")

# Clean and subset to CS snps with PIP>0.8
# For each row, filter cs_snps by corresponding cs_pips > 0.8, then combine all rsIDs

rsids_high_pip <- unique(unlist(
  lapply(1:nrow(dat), function(i) {
    snps <- dat$cs_snps[[i]]
    pips <- dat$cs_pips[[i]]
    snps[pips > 0.8]
  })
))

# Get eQTL and GWAS rsIDs
rsids_high_pip <- c(pvar$V2[which(pvar$id1 %in% rsids_high_pip)], pvar$V2[which(pvar$id2 %in% rsids_high_pip)])

# # Now divide into 10 bins based on phe_id (gene-level ID)
# gene_list <- sort(unique(dat$phe_id))
# n_genes <- length(gene_list)

# # Assign each gene to a bin
# bins <- cut(seq_along(gene_list), breaks = 10, labels = FALSE)

# # Ensure chunk is valid
# if (chunk < 1 || chunk > 10) stop("Invalid chunk number: must be between 1 and 100")

# # Select genes in the chosen bin
# genes_in_chunk <- gene_list[bins == chunk]

# Subset dat to just those genes
# dat_chunk <- dat[dat$phe_id %in% genes_in_chunk, ]

# pvar <- data.frame(fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr,".bim")))
# pvar <- pvar[pvar$V4 %in% dat_chunk$var_from,]
# pvar$chr_pos <- paste0("chr",pvar$V1,":",pvar$V2)

# Get eQTL and GWAS rsIDs
eqtl_rsids <- unique(na.omit(rsids_high_pip))
gwas_rsids <- unique(na.omit(dat_bc$RSID))
all_rsids <- unique(c(eqtl_rsids, gwas_rsids))

# Define inputs
plink_prefix <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr)  # without file extensions
out_prefix <- paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/gwas_loci/TCGA_finemap/ld_results_chr",chr,"_",assembly,"_",level,"_k",k,"_L",L)
snp_file <- tempfile(pattern = "snplist_", fileext = ".txt")
fwrite(data.table(SNP = all_rsids), snp_file, col.names = FALSE, quote = FALSE)

# Construct PLINK2 command
plink_cmd <- paste(
  "plink2",
  "--bfile", plink_prefix,
  "--r2-unphased",
  "--extract", snp_file,
  "--ld-window", 99999,
  "--ld-window-kb", 1000,
  "--ld-window-r2", 0.01,
  "--out", out_prefix
)

# Run the command
system(plink_cmd)

# Example path to your .vcor file

# Read it in
vcor_dat <- fread(paste0(out_prefix,".vcor"))

# LD between eQTLs and GWAS SNPs
ld_eqtl_gwas1 <- vcor_dat[
  (ID_A %in% eqtl_rsids & ID_B %in% gwas_rsids)
]
tmp1 <- ld_eqtl_gwas1
colnames(tmp1) <-c ("eqtl_chr","eqtl_pos","eqtl_rsid","gwas_chr","gwas_pos","gwas_rsid","r2")

ld_eqtl_gwas2 <- vcor_dat[
  (ID_B %in% eqtl_rsids & ID_A %in% gwas_rsids)
]
tmp2 <- ld_eqtl_gwas2[,c(4,5,6,1,2,3,7)]
colnames(tmp2) <-c ("eqtl_chr","eqtl_pos","eqtl_rsid","gwas_chr","gwas_pos","gwas_rsid","r2")
out <- rbind(tmp1,tmp2)


file.remove(snp_file)

write.table(out,file=paste0(out_prefix,".olap.PIP0.8.CS"),row.names=F,col.names=T,sep="\t",quote=F)


##############################
# now do all CS SNPs regardless of PIP

rsids_no_pip <- unique(unlist(
  lapply(1:nrow(dat), function(i) {
    snps <- dat$cs_snps[[i]]
  })
))

rsids_no_pip <- c(pvar$V2[which(pvar$id1 %in% rsids_no_pip)], pvar$V2[which(pvar$id2 %in% rsids_no_pip)])

# Get eQTL and GWAS rsIDs
eqtl_rsids <- unique(na.omit(rsids_no_pip))
gwas_rsids <- unique(na.omit(dat_bc$RSID))
all_rsids <- unique(c(eqtl_rsids, gwas_rsids))

# Define inputs
plink_prefix <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr)  # without file extensions
out_prefix <- paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/gwas_loci/TCGA_finemap/ld_results_chr",chr,"_",assembly,"_",level,"_k",k,"_L",L)
snp_file <- tempfile(pattern = "snplist_", fileext = ".txt")
fwrite(data.table(SNP = all_rsids), snp_file, col.names = FALSE, quote = FALSE)

# Construct PLINK2 command
plink_cmd <- paste(
  "plink2",
  "--bfile", plink_prefix,
  "--r2-unphased",
  "--extract", snp_file,
  "--ld-window", 99999,
  "--ld-window-kb", 1000,
  "--ld-window-r2", 0.01,
  "--out", out_prefix
)

# Run the command
system(plink_cmd)

# Example path to your .vcor file

# Read it in
vcor_dat <- fread(paste0(out_prefix,".vcor"))

# LD between eQTLs and GWAS SNPs
ld_eqtl_gwas1 <- vcor_dat[
  (ID_A %in% eqtl_rsids & ID_B %in% gwas_rsids)
]
tmp1 <- ld_eqtl_gwas1
colnames(tmp1) <-c ("eqtl_chr","eqtl_pos","eqtl_rsid","gwas_chr","gwas_pos","gwas_rsid","r2")

ld_eqtl_gwas2 <- vcor_dat[
  (ID_B %in% eqtl_rsids & ID_A %in% gwas_rsids)
]
tmp2 <- ld_eqtl_gwas2[,c(4,5,6,1,2,3,7)]
colnames(tmp2) <-c ("eqtl_chr","eqtl_pos","eqtl_rsid","gwas_chr","gwas_pos","gwas_rsid","r2")
out <- rbind(tmp1,tmp2)


file.remove(snp_file)

write.table(out,file=paste0(out_prefix,".olap.PIPNA.CS"),row.names=F,col.names=T,sep="\t",quote=F)




