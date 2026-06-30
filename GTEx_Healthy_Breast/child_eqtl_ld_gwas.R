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

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
assembly <- as.character(args[1])
level <- as.character(args[2])
k <- as.numeric(args[3])
eqtl_dir <- as.character(args[4]) 
chr <- as.numeric(args[5])
chunk <- as.numeric(args[6])

# assembly="gencodev45"
# level="gene"
# k=50
# eqtl_dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools"
# chr=1
# chunk=1

####################################################################################
# begin code
####################################################################################

# Load and subset GWAS loci
dat_bc <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/gwas_loci/indep_BC_gwas_loci_332.txt"))
dat_bc <- dat_bc[which(dat_bc$CHR == chr), ]

if (nrow(dat_bc) == 0) {
  stop("No GWAS loci on this chromosome — exiting job.")
}

# Load conditional eQTLs
file <- paste0(eqtl_dir, "/GTEx_cond_", assembly, "_ODS1_TMM1_VST1_lstpm0_doBC0_k", k, "_HC_genpc5/", level, "_qtls_cond_nperm1000_normalized_full.txt.gz")
dat <- fread(file, header = FALSE)

# Load and apply header
header <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/conditional_header_0.txt", header = T)
names(dat) <- names(header)

# Clean and subset to significant QTLs on chromosome
dat <- dat[dat$bwd_sig == 1, ]
dat$assembly <- assembly
dat$SNP_clean <- paste(dat$var_chr, dat$var_from, sep = ":")
dat <- dat[which(dat$phe_chr == paste0("chr", chr)), ]

# Now divide into 10 bins based on phe_id (gene-level ID)
gene_list <- sort(unique(dat$phe_id))
n_genes <- length(gene_list)

# Assign each gene to a bin
bins <- cut(seq_along(gene_list), breaks = 10, labels = FALSE)

# Ensure chunk is valid
if (chunk < 1 || chunk > 10) stop("Invalid chunk number: must be between 1 and 100")

# Select genes in the chosen bin
genes_in_chunk <- gene_list[bins == chunk]

# Subset dat to just those genes
dat_chunk <- dat[dat$phe_id %in% genes_in_chunk, ]
# dat_chunk$var_id2 <- paste(dat_chunk$SNP_clean,dat_chunk$)

# pvar <- data.frame(fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr,".bim")))
# pvar <- pvar[pvar$V4 %in% dat_chunk$var_from,]
# # pvar$chr_pos <- paste0("chr",pvar$V1,":",pvar$V2)
# pvar$id1 <- paste(pvar$V1,pvar$V4,pvar$V5,pvar$V6,sep=":")
# pvar$id2 <- paste(pvar$V1,pvar$V4,pvar$V6,pvar$V5,sep=":")

# Get eQTL and GWAS rsIDs
eqtl_rsids <- unique(dat_chunk$var_id)
gwas_rsids <- unique(na.omit(dat_bc$RSID))
all_rsids <- unique(c(eqtl_rsids, gwas_rsids))

# Define inputs
plink_prefix <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/EUR/chr",chr)  # without file extensions
out_prefix <- paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/gwas_loci/GTEx_eqtl/ld_results_chr",chr,"_chunk",chunk,"_",assembly,"_",level,"_k",k)
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

write.table(out,file=paste0(out_prefix,".olap"),row.names=F,col.names=T,sep="\t",quote=F)



