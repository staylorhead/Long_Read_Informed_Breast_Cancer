#!/usr/bin/env Rscript

##################################################################################
# change library to local
##################################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

####################################################################################
# load dependencies
####################################################################################
library(MungeSumstats)
library(data.table)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
ref_genome <- as.character(args[1])
target_genome <- as.character(args[2])
pheno <- as.character(args[3])

####################################################################################
# begin
####################################################################################

ms_dat <- fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,".tsv.gz"))
ms_dat <- data.frame(ms_dat)

lifted <- liftover(
  sumstats_dt=ms_dat,
  convert_ref_genome=target_genome,
  ref_genome=ref_genome,
  chain_source = "ensembl",
  imputation_ind = FALSE,
  chrom_col = "CHR",
  start_col = "BP",
  end_col = "BP",
  as_granges = FALSE,
  style = "NCBI",
  verbose = TRUE
)

fwrite(lifted, file=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",target_genome,"_",pheno,".tsv.gz"))

rm(lifted)

ms_dat <- fread(paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed.tsv.gz"))
ms_dat <- data.frame(ms_dat)

lifted <- liftover(
  sumstats_dt=ms_dat,
  convert_ref_genome=target_genome,
  ref_genome=ref_genome,
  chain_source = "ensembl",
  imputation_ind = FALSE,
  chrom_col = "CHR",
  start_col = "BP",
  end_col = "BP",
  as_granges = FALSE,
  style = "NCBI",
  verbose = TRUE
)

fwrite(lifted, file=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",target_genome,"_",pheno,"_ambig_strand_removed.tsv.gz"))


####################################################################################
# check memory usage
####################################################################################
print(gc())
