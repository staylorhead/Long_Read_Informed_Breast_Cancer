#!/usr/bin/env Rscript

########################################################################
# change library to local
########################################################################
myPaths <- .libPaths()
myPaths <- c("/rsrch5/home/epi/sthead/R/x86_64-pc-linux-gnu-library/4.3",myPaths)
.libPaths(myPaths)

########################################################################
# load dependencies
########################################################################
library(data.table)
library(dplyr)

# clear out environment
rm(list = ls())

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
assembly <- as.character(args[1]) # gencodev45 or veigaTum
do_ods <- as.numeric(args[2]) # 0 or 1
lstpm <- as.numeric(args[3]) # 0 or 1
do_vst <- as.numeric(args[4]) # 0 or 1
do_tmm <- as.numeric(args[5]) # 0 or 1
do_bc <- as.numeric(args[6])# 0 or 1
k <- as.numeric(args[7])
method <- as.character(args[8]) # "HC" or "PC"
gen_pcs_file <- as.character(args[9])
num_gen_pcs <- as.numeric(args[10])

# assembly="gencodev45"
# do_ods=1
# lstpm=0
# do_vst=1
# do_tmm=1
# do_bc=0
# k=15
# method="HC"
# gen_pcs_file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx_838_pca.eigenvec"
# num_gen_pcs=5

########################################################################
# helper functions
########################################################################
hidden_convariate_linear <- function(F,Y,k,lambda,lambda2,lambda3,iter) {
  ## Use Example
  # hcp = hidden_convariate_linear(standardize(datSeq), standardize(t(datExpr)),k=10,iter = 100)
  #
  # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
  # input:
  #      F: a matrix nxd of known covariates, where n is the number of
  #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
  #      Y: a matrix of nxg of expression data (must be standardized (columns
  #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
  #      k: number of inferred hidden components (k is an integer)
  #      lambda, lambda2, lambda3 are model parameters
  #      (optional) iter: number of iterations (default = 100);
  #
  #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
  #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
  #      typically, if lambda>5, then hidden factors match the known covariates closely.
  #
  # objective:
  #
  # this function solves the following problem:
  # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
  #
  # output:
  #      Z: matrix of hidden components, dimensionality: nxk
  #      B: matrix of effects of hidden components, dimensionality: kxg
  #      o: value of objective function on consecutive iterations.
  #
  # to use the residual data: Residual = Y - Z*B
  library(MASS)
  library(pracma)
  
  tol = 1e-6;
  
  U = matrix(0, nrow=dim(F)[2],k)
  Z = matrix(0, nrow=dim(F)[1],k)
  B = matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
  F = as.matrix(F)
  
  n1 = dim(F)[1]
  d1 = dim(F)[2]
  
  n2 = dim(Y)[1]
  d2 = dim(Y)[2]
  
  if(n1!=n2)    stop("number of rows in F and Y must agree")
  
  if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
    stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
  }
  
  o = vector(length=iter)
  
  for (ii in 1:iter) {
    print(ii)
    o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + 
      (sum(B^2))*lambda2 + lambda3*(sum(U^2));
    Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
    B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
    U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)
    
    if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
  }
  
  error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
  error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
  error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);
  
  dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
  db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
  du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;
  
  
  dataout = list(Z = Z, B = B, U = U)
  return(dataout)
}

standardize<- function(X){
  X = as.matrix(X)
  # n = dim(X)[1]
  # p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale = F)
  # X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))
  
  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  return (X)
}

########################################################################
# begin HC estimation or expression PC loading
########################################################################

set.seed(555)

prefix <- paste0(assembly,"_ODS",do_ods,"_TMM",do_tmm,"_VST",do_vst,"_lstpm",lstpm,"_doBC",do_bc)

gene_bed = fread(paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_gene_exp_',prefix,'.bed.gz'))
gene_cbt = as.matrix(gene_bed[,7:ncol(gene_bed)])

tx_bed = fread(paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_tx_exp_',prefix,'.bed.gz'))
tx_cbt = as.matrix(tx_bed[,7:ncol(tx_bed)])

sum(colnames(gene_cbt)==colnames(tx_cbt))==ncol(gene_cbt) # good
tpm_ids <- colnames(gene_cbt)

# read in prior data
setwd("/rsrch5/scratch/epi/sthead/GTEx_v8/picard_multiqc/multiqc_data")

if(method=="HC"){
  dat_rna <- fread("multiqc_picard_RnaSeqMetrics.txt")
  dat_summ <- fread("multiqc_general_stats.txt")
  dat_gen <- fread("multiqc_picard_AlignmentSummaryMetrics.txt")

  sample_dat <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"))
  sample_dat$gtex_id <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", sample_dat$sample)
  rm <- which(sample_dat$SMAFRZE=="EXCLUDE")
  sample_dat <- sample_dat[-rm,]
  rownames(sample_dat) <- sample_dat$gtex_id
  keep_ids_full <- sample_dat$sample

  sum(dat_rna$Sample==dat_summ$Sample) # good
  sum(dat_rna$Sample==dat_gen$Sample) # good

  dat_prior <- cbind(dat_rna, dat_summ[,-1],dat_gen[,-1])
  dat_prior <- data.frame(dat_prior)
  dat_prior <- dat_prior[dat_prior$Sample %in% keep_ids_full,]
  dat_prior$Sample <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", dat_prior$Sample)
  dat_prior <- dat_prior[dat_prior$Sample %in% tpm_ids,]
  
  prior_barcodes <- dat_prior$Sample
  match_in_prior <- sapply(tpm_ids,FUN=function(x){
    which(prior_barcodes==x)
  })
  dat_prior <- dat_prior[match_in_prior,]

  if(sum(dat_prior$Sample==tpm_ids)==length(tpm_ids)){
        # remove non-numeric columns
      numeric_cols <- sapply(dat_prior, is.numeric)
      keep <- which(numeric_cols==T)
      dat_prior_sub <-  dat_prior[, keep]

      # remove columns with no variation
      dat_prior_sub <- dat_prior_sub[, sapply(dat_prior_sub, function(x) length(unique(x)) > 1), drop=FALSE]

      hcp_gene <- hidden_convariate_linear(standardize(dat_prior_sub), 
                                      standardize(t(gene_cbt)), 
                                      lambda=5,lambda2=1, lambda3=1, 
                                      k=k, iter=100)

      hcp_tx <- hidden_convariate_linear(standardize(dat_prior_sub), 
                                      standardize(t(tx_cbt)), 
                                      lambda=5,lambda2=1, lambda3=1, 
                                      k=k, iter=100)
    }else{
      stop("error in sample matching")
    }
  }else if(method=="PC"){

    # FIX THIS LATER. JUST RUN WITH METHOD=HC
  load(paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/pca/expr_pca_dat_",prefix,".RData"))
  pc_gene <- data.frame(pca_gene$x)
  pc_tx <- data.frame(pca_tx$x)

  sum(colnames(gene_cbt)==rownames(pc_gene)) # good
  sum(colnames(tx_cbt)==rownames(pc_tx)) # good
  pc_gene <- pc_gene[,1:k]
  pc_tx <- pc_tx[,1:k]

}

### aggregate covariates

individual <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
individual$SEX <- as.character(individual$SEX)
individual$sex_out = factor(individual$SEX,levels=c("2","1"),labels=c("female","male"))
individual$age_out <- sapply(strsplit(individual$AGE, "-"), function(x) {
  mean(as.numeric(x))
})
covariates <- individual
rownames(covariates) <- covariates$SUBJID
covariates <- covariates[,c("SUBJID","sex_out","age_out")]
colnames(covariates) <- c("SUBJID","sex","age")
#covariates$age2 = covariates$age^2

if(method=="HC"){
  hcp_ids <- rownames(hcp_gene$Z)
  covariates <- covariates[hcp_ids,]
  sum(covariates$SUBJID==hcp_ids) # good

  covariates_gene = cbind(covariates,hcp_gene$Z)
  covariates_tx = cbind(covariates, hcp_tx$Z)

  }else if(method=="PC"){
    ## will need to fix this eventually
  sum(covariates$individualID==rownames(pc_gene)) # good
  covariates_gene = cbind(covariates,pc_gene)
  covariates_tx = cbind(covariates,pc_tx)

}

colnames(covariates_gene)[-1:-3] = paste0('hcp',1:k)
colnames(covariates_tx)[-1:-3] = paste0('hcp',1:k)

gen_pc = data.frame(fread(gen_pcs_file))
gen_pc$V1 <- gsub("\\.", "-", gen_pc$V1)
gen_pc <- gen_pc[,-2]
gen_pc <- gen_pc[,1:(num_gen_pcs+1)] # keep first X pcs
colnames(gen_pc)[-1] <- paste0("PC",1:num_gen_pcs)

length(unique(gen_pc$V1)) # good

covariates_gene <- merge(x=covariates_gene,y=gen_pc,by.x="SUBJID",by.y="V1",all.x=T,all.y=F)
covariates_tx <- merge(x=covariates_tx,y=gen_pc,by.x="SUBJID",by.y="V1",all.x=T,all.y=F)

dim(covariates_gene)
dim(covariates_tx)

# restrict samples to only those with good genotypes
covariates_gene <- covariates_gene[!is.na(covariates_gene$PC1),]
covariates_tx <- covariates_tx[!is.na(covariates_tx$PC1),]
dim(covariates_gene)
dim(covariates_tx)

covariates_gene$sex <- as.numeric(covariates_gene$sex)-1
covariates_tx$sex <- as.numeric(covariates_tx$sex)-1

colnames(covariates_gene) <- c("individualID","sex","age",paste0("hcp",1:k),paste0("PC",1:num_gen_pcs))
cov_gene <- covariates_gene
colnames(covariates_tx) <- c("individualID","sex","age",paste0("hcp",1:k),paste0("PC",1:num_gen_pcs))
cov_tx <- covariates_tx

tcov_gene <- t(cov_gene)
row.names(tcov_gene)[1] <- "id"

tcov_tx <- t(cov_tx)
row.names(tcov_tx)[1] <- "id"

########################################################################
# write out files
########################################################################

write.table(tcov_gene,file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_cov_qtltools_gene_",prefix,"_k",k,"_",method,"_genpc",num_gen_pcs,".txt"),
  sep="\t",row.names=T,col.names=F,quote=F)

write.table(tcov_tx,file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_cov_qtltools_tx_",prefix,"_k",k,"_",method,"_genpc",num_gen_pcs,".txt"),
  sep="\t",row.names=T,col.names=F,quote=F)

