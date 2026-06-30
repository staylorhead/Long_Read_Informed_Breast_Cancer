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
library(readr)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
ref_genome <- as.character(args[1])
pheno <- as.character(args[2])

####################################################################################
# begin
####################################################################################

## OVERALL BC

if(pheno=="overallBC"){

dat_big <- read_table("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics-003.txt")
dat <- dat_big[,c("var_name","EAFcontrols.iCOGs","Effect.iCOGs","r2.iCOGs","r2.Onco","chr.iCOGs","Position.iCOGs","Effect.Meta","Baseline.Meta","Beta.meta",
  "var.meta","sdE.meta","chi.meta","p.meta")]
dat <- na.omit(dat)
sum(dat$Effect.iCOGs==dat$Effect.Meta,na.rm=T) # ok
# all Onco and iCOGs imputation R2 greater than 0.3
dat$minR2 <- pmin(dat$r2.iCOGs, dat$r2.Onco, na.rm = TRUE) 
dat <- dat[,c("var_name","EAFcontrols.iCOGs","minR2","chr.iCOGs","Position.iCOGs","Effect.Meta","Baseline.Meta","Beta.meta",
  "var.meta","sdE.meta","chi.meta","p.meta")]
colnames(dat) <- c("SNP","ALTERNATIVE_AF","INFO","CHR","BP","EFFECT.ALLELE","REFERENCE.ALLELE","BETA","VAR","SE","CHI","P")
dat <- data.table(dat)
dat$BP <- as.integer(dat$BP)
dat$EFFECT.ALLELE <- toupper(dat$EFFECT.ALLELE)
dat$REFERENCE.ALLELE <- toupper(dat$REFERENCE.ALLELE)

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,".tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5
# )

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed.tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5,
#   strand_ambig_filter=TRUE
# )

ms_dat <- format_sumstats(
  path=dat,
  ref_genome=ref_genome,
  save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed_no_biallelic_filt.tsv.gz"),
  INFO_filter = 0.3,
  nThread = 5,
  strand_ambig_filter=TRUE,
  bi_allelic_filter=FALSE
)

}else if(pheno %in% c("lumA","lumB","lumBHER2neg","HER2enriched")){


####### SUBTYPE

dat_big <- fread("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics-002.txt")

if(pheno=="lumA"){

  dat <- dat_big[,c("var_name","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs","Effect.iCOGs","r2.iCOGs","r2.Onco","Effect.Meta","Baseline.Meta",
    "Luminal_A_log_or_meta","Luminal_A_se_meta")]
  colnames(dat)[10:11] <- c("BETA","SE")

  }else if(pheno=="lumB"){

    dat <- dat_big[,c("var_name","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs","Effect.iCOGs","r2.iCOGs","r2.Onco","Effect.Meta","Baseline.Meta",
    "Luminal_B_log_or_meta","Luminal_B_se_meta")]
    colnames(dat)[10:11] <- c("BETA","SE")

    }else if(pheno=="lumBHER2neg"){

      dat <- dat_big[,c("var_name","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs","Effect.iCOGs","r2.iCOGs","r2.Onco","Effect.Meta","Baseline.Meta",
      "Luminal_B_HER2Neg_log_or_meta","Luminal_B_HER2Neg_se_meta")]
      colnames(dat)[10:11] <- c("BETA","SE")

      }else if(pheno=="HER2enriched"){
        dat <- dat_big[,c("var_name","chr.iCOGs","Position.iCOGs","EAFcontrols.iCOGs","Effect.iCOGs","r2.iCOGs","r2.Onco","Effect.Meta","Baseline.Meta",
        "HER2_Enriched_log_or_meta","HER2_Enriched_se_meta")]
        colnames(dat)[10:11] <- c("BETA","SE")

}

dat <- na.omit(dat)
dat$minR2 <- pmin(dat$r2.iCOGs, dat$r2.Onco, na.rm = TRUE) # all over 0.3

dat <- dat[,c("var_name","EAFcontrols.iCOGs","chr.iCOGs","Position.iCOGs","Effect.Meta","Baseline.Meta","BETA","SE","minR2")]
colnames(dat) <- c("SNP","ALTERNATIVE_AF","CHR","BP","EFFECT.ALLELE","REFERENCE.ALLELE","BETA","SE","INFO")
dat <- data.table(dat)
dat$BP <- as.integer(dat$BP)
dat$EFFECT.ALLELE <- toupper(dat$EFFECT.ALLELE)
dat$REFERENCE.ALLELE <- toupper(dat$REFERENCE.ALLELE)

z = dat$BETA / dat$SE
dat$P = 2 * pnorm(-abs(z))

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,".tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5
# )

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed.tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5,
#   strand_ambig_filter=TRUE
# )  

ms_dat <- format_sumstats(
  path=dat,
  ref_genome=ref_genome,
  save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed_no_biallelic_filt.tsv.gz"),
  INFO_filter = 0.3,
  nThread = 5,
  strand_ambig_filter=TRUE,
  bi_allelic_filter=FALSE
)  




}else if(pheno=="TN"){

####### CIMBA TRIPLE NEGATIVE

dat_big <- fread("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt")

dat <- na.omit(dat_big)

dat <- dat[,c("MarkerName","eff_freq","CHR","position","Allele2","Allele1","Effect","StdErr","P-value")]
colnames(dat) <- c("SNP","ALTERNATIVE_AF","CHR","BP","EFFECT.ALLELE","REFERENCE.ALLELE","BETA","SE","P")
dat <- data.table(dat)
dat$BP <- as.integer(dat$BP)
dat$EFFECT.ALLELE <- toupper(dat$EFFECT.ALLELE)
dat$REFERENCE.ALLELE <- toupper(dat$REFERENCE.ALLELE)

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,".tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5
# )

# ms_dat <- format_sumstats(
#   path=dat,
#   ref_genome=ref_genome,
#   save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed.tsv.gz"),
#   INFO_filter = 0.3,
#   nThread = 5,
#   strand_ambig_filter=TRUE
# )  

ms_dat <- format_sumstats(
  path=dat,
  ref_genome=ref_genome,
  save_path=paste0("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_",ref_genome,"_",pheno,"_ambig_strand_removed_no_biallelic_filt.tsv.gz"),
  INFO_filter = 0.3,
  nThread = 5,
  strand_ambig_filter=TRUE,
  bi_allelic_filter=FALSE
)  

}


####################################################################################
# check memory usage
####################################################################################
print(gc())
