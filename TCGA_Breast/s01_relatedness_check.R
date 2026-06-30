# load dependencies
library(SummarizedExperiment)
library(data.table)
library(readxl)
library(dplyr)
library(tximeta)

# load tximeta object gencode v38 for RNA-seq sample IDs
load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/se_gencodev38.RData")
sample_ids_se <- colnames(assay(se))
length(sample_ids_se)
# 1231
rm(se) 

# extract sample information from RNA-seq barcodes
sample_dat <- data.frame(barcode=sample_ids_se)
sample_dat$TSS <- unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[2]]))
sample_dat$participant <- unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[3]]))
sample_dat$rna_sample_type <- substr(unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[4]])),1,2)
sample_dat$participant <- substr(sample_dat$barcode,1,12)
sample_dat$rna_sample_type <- factor(sample_dat$rna_sample_type,levels=c("01","06","11"),
  labels=c("Primary Solid Tumor","Metastatic","Solid Tissue Normal"))

# read in sample IDs from genotype data
geno_dat <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC/allchr.smiss")
nrow(geno_dat)
# 1150
colnames(geno_dat)[1] <- "IID"
geno_dat$participant <- substr(geno_dat$IID,1,12)
geno_dat$sample_type <- substr(unlist(lapply(strsplit(geno_dat$IID,"-"), function(l) l[[4]])),1,2)
geno_dat$sample_type <- factor(geno_dat$sample_type,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))

# read in subtype data
subtype_dat <- read_excel("/rsrch5/home/epi/sthead/TCGA/covariates/subtype_data.xlsx")

# read in genotype relatedness table (no imputation R2 filtering)
king_dat <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC/allchr.kin0")
colnames(king_dat)[1] <- "IID"
king_dat$participant1 <- substr(king_dat$IID,1,12)
king_dat$participant2 <- substr(king_dat$IID2,1,12)

king_dat$sample_type1 <- substr(unlist(lapply(strsplit(king_dat$IID,"-"), function(l) l[[4]])),1,2)
king_dat$sample_type1 <- factor(king_dat$sample_type1,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))
king_dat$sample_type2 <- substr(unlist(lapply(strsplit(king_dat$IID2,"-"), function(l) l[[4]])),1,2)
king_dat$sample_type2 <- factor(king_dat$sample_type2,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))

related <- king_dat[king_dat$KINSHIP>0.0884,]
dim(related)
# 54

# double check that the same pairs aren't repeated in this table
related_unique <- related %>%
  mutate(pair_id = pmin(IID, IID2), pair_id2 = pmax(IID, IID2)) %>%
  distinct(pair_id, pair_id2, .keep_all = TRUE) %>%
  select(-pair_id, -pair_id2)  # 
dim(related_unique)
# 54
choose(1150,2)==nrow(king_dat)
# TRUE
# no repeat entries in related dataframe (same pair of IDs)

sum(related$participant1==related$participant2)
# 54
# all related pairs are different samples from the same individual

# initiate vector of unrelated samples to proceed with
good_ids <- geno_dat$IID[!geno_dat$IID %in% related$IID & !geno_dat$IID %in% related$IID2]
length(good_ids)
# [1] 1042
1042+54*2
# [1] 1150

table(paste0(related$sample_type1,related$sample_type2))
# Blood Derived NormalBlood Derived Normal
#                                        2
#  Solid Tissue NormalBlood Derived Normal
#                                       52
# either two blood draws from the same person or tissue + blood draw for same person

twobloods <- related[related$sample_type1=="Blood Derived Normal" & 
                       related$sample_type2=="Blood Derived Normal",]
dim(twobloods)
# [1]  2 10
blood_tiss <- related[related$sample_type1=="Solid Tissue Normal" & 
                        related$sample_type2=="Blood Derived Normal",]
dim(blood_tiss)
# [1] 52 10

# for those with two blood samples:
# preferentially keep blood sample with lower missingness
for(i in 1:nrow(twobloods)){
	samp1 <- twobloods$IID[i]
	samp2 <- twobloods$IID2[i]
	miss1 <- geno_dat$F_MISS[geno_dat$IID==samp1]
	miss2 <- geno_dat$F_MISS[geno_dat$IID==samp2]
	keep <- which.min(c(miss1,miss2))
	good_ids <- c(good_ids, c(samp1,samp2)[keep])
	rm(keep)
}

# for those with one tissue and one blood sample:
# keep the blood sample
for(i in 1:nrow(blood_tiss)){
	samp1 <- blood_tiss$IID[i]
	samp2 <- blood_tiss$IID2[i]
	keep <- which(c(blood_tiss$sample_type1[i],blood_tiss$sample_type2[i])=="Blood Derived Normal")
	good_ids <- c(good_ids, c(samp1,samp2)[keep])
	rm(keep)
}
length(good_ids)
#[1] 1096

tmp <- king_dat[king_dat$IID %in% good_ids & king_dat$IID2 %in% good_ids,]
summary(tmp$KINSHIP)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.551223 -0.173221 -0.018826 -0.084041 -0.001473  0.051155

# save list of unrelated blood samples
out <- geno_dat[geno_dat$IID %in% good_ids,]
dim(out)
# 1096 unrelated samples

table(out$sample_type)
# Blood Derived Normal  Solid Tissue Normal
#                 1009                   87

length(unique(out$participant))
# 1096
# one sample per person and minimal relatedness
unrelated_samples_no_impR2_filt <- good_ids

write.table(out$IID,file="/rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt",
	sep="\t",col.names=F,row.names=F,quote=F)


## now, with imputation R2 filter


# load tximeta object gencode v38 for RNA-seq sample IDs
load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/se_gencodev38.RData")
sample_ids_se <- colnames(assay(se))
length(sample_ids_se)
# 1231
rm(se) 

# extract sample information from RNA-seq barcodes
sample_dat <- data.frame(barcode=sample_ids_se)
sample_dat$TSS <- unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[2]]))
sample_dat$participant <- unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[3]]))
sample_dat$rna_sample_type <- substr(unlist(lapply(strsplit(sample_dat$barcode,"-"), function(l) l[[4]])),1,2)
sample_dat$participant <- substr(sample_dat$barcode,1,12)
sample_dat$rna_sample_type <- factor(sample_dat$rna_sample_type,levels=c("01","06","11"),
  labels=c("Primary Solid Tumor","Metastatic","Solid Tissue Normal"))

# read in sample IDs from genotype data
geno_dat <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03/allchr.smiss")
nrow(geno_dat)
# 1150
colnames(geno_dat)[1] <- "IID"
geno_dat$participant <- substr(geno_dat$IID,1,12)
geno_dat$sample_type <- substr(unlist(lapply(strsplit(geno_dat$IID,"-"), function(l) l[[4]])),1,2)
geno_dat$sample_type <- factor(geno_dat$sample_type,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))

# read in subtype data
subtype_dat <- read_excel("/rsrch5/home/epi/sthead/TCGA/covariates/subtype_data.xlsx")

# read in genotype relatedness table (no imputation R2 filtering)
king_dat <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03/allchr.kin0")
colnames(king_dat)[1] <- "IID"
king_dat$participant1 <- substr(king_dat$IID,1,12)
king_dat$participant2 <- substr(king_dat$IID2,1,12)

king_dat$sample_type1 <- substr(unlist(lapply(strsplit(king_dat$IID,"-"), function(l) l[[4]])),1,2)
king_dat$sample_type1 <- factor(king_dat$sample_type1,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))
king_dat$sample_type2 <- substr(unlist(lapply(strsplit(king_dat$IID2,"-"), function(l) l[[4]])),1,2)
king_dat$sample_type2 <- factor(king_dat$sample_type2,levels=c("10","11"),
  labels=c("Blood Derived Normal","Solid Tissue Normal"))

related <- king_dat[king_dat$KINSHIP>0.0884,]
dim(related)
# 54

# double check that the same pairs aren't repeated in this table
related_unique <- related %>%
  mutate(pair_id = pmin(IID, IID2), pair_id2 = pmax(IID, IID2)) %>%
  distinct(pair_id, pair_id2, .keep_all = TRUE) %>%
  select(-pair_id, -pair_id2)  # 
dim(related_unique)
# 54
choose(1150,2)==nrow(king_dat)
# TRUE
# no repeat entries in related dataframe (same pair of IDs)

sum(related$participant1==related$participant2)
# 54
# all related pairs are different samples from the same individual

# initiate vector of unrelated samples to proceed with
good_ids <- geno_dat$IID[!geno_dat$IID %in% related$IID & !geno_dat$IID %in% related$IID2]
length(good_ids)
# [1] 1042
1042+54*2
# [1] 1150

table(paste0(related$sample_type1,related$sample_type2))
# Blood Derived NormalBlood Derived Normal
#                                        2
#  Solid Tissue NormalBlood Derived Normal
#                                       52
# either two blood draws from the same person or tissue + blood draw for same person

twobloods <- related[related$sample_type1=="Blood Derived Normal" & 
                       related$sample_type2=="Blood Derived Normal",]
dim(twobloods)
# [1]  2 10
blood_tiss <- related[related$sample_type1=="Solid Tissue Normal" & 
                        related$sample_type2=="Blood Derived Normal",]
dim(blood_tiss)
# [1] 52 10

# for those with two blood samples:
# preferentially keep blood sample with lower missingness
for(i in 1:nrow(twobloods)){
	samp1 <- twobloods$IID[i]
	samp2 <- twobloods$IID2[i]
	miss1 <- geno_dat$F_MISS[geno_dat$IID==samp1]
	miss2 <- geno_dat$F_MISS[geno_dat$IID==samp2]
	keep <- which.min(c(miss1,miss2))
	good_ids <- c(good_ids, c(samp1,samp2)[keep])
	rm(keep)
}

# for those with one tissue and one blood sample:
# keep the blood sample
for(i in 1:nrow(blood_tiss)){
	samp1 <- blood_tiss$IID[i]
	samp2 <- blood_tiss$IID2[i]
	keep <- which(c(blood_tiss$sample_type1[i],blood_tiss$sample_type2[i])=="Blood Derived Normal")
	good_ids <- c(good_ids, c(samp1,samp2)[keep])
	rm(keep)
}
length(good_ids)
#[1] 1096

tmp <- king_dat[king_dat$IID %in% good_ids & king_dat$IID2 %in% good_ids,]
summary(tmp$KINSHIP)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.549984 -0.172632 -0.018832 -0.083793 -0.001471  0.051195

# save list of unrelated blood samples
out <- geno_dat[geno_dat$IID %in% good_ids,]
dim(out)
# 1096 unrelated samples

table(out$sample_type)
# Blood Derived Normal  Solid Tissue Normal
#                 1009                   87

length(unique(out$participant))
# 1096
# one sample per person and minimal relatedness
sum(good_ids %in% unrelated_samples_no_impR2_filt)
# 1096, same results as with no imputation R2 filtering, good
