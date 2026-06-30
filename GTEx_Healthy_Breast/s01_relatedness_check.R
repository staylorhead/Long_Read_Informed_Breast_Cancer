# load dependencies
library(SummarizedExperiment)
library(data.table)
library(readxl)
library(dplyr)
library(tximeta)

# read in genotype relatedness table (no imputation R2 filtering)
king_dat <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/geno/GTEx_MAF0.01_passQC/allchr.kin0")
#colnames(king_dat)[1] <- "IID"
king_dat$participant1 <- king_dat$IID1
king_dat$participant2 <- king_dat$IID2

related <- king_dat[king_dat$KINSHIP>0.0884,]
dim(related)
# 0

summary(king_dat$KINSHIP)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.252538 -0.039604 -0.005971 -0.031446 -0.001331  0.067101