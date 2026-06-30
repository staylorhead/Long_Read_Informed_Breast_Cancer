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
library(tximeta)
library(fishpond)
library(SummarizedExperiment)
library(GenomicFeatures)
library(VariantAnnotation)
library(rtracklayer)
library(Biostrings)
library(stringr)
library(VGAM)
library(GenomicRanges)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(dplyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(sva)

# clear out environment
rm(list = ls())

########################################################################
# parse arguments
########################################################################
args <- commandArgs(trailingOnly = TRUE)
assembly <- as.character(args[1]) # gencodev45 or veigaNorm or comb
do_ods <- as.numeric(args[2]) # 0 or 1
lstpm <- as.numeric(args[3]) # 0 or 1
do_vst <- as.numeric(args[4]) # 0 or 1
do_tmm <- as.numeric(args[5]) # 0 or 1
do_bc <- as.numeric(args[6])# 0 or 1

########################################################################
# helper functions
########################################################################
makeCountsFromAbundance <- function(countsMat, abundanceMat, lengthMat,
                                    countsFromAbundance=c("scaledTPM","lengthScaledTPM")) {
  #countsFromAbundance <- match.arg(countsFromAbundance)
  sparse <- is(countsMat, "dgCMatrix")
  colsumfun <- if (sparse) Matrix::colSums else colSums
  countsSum <- colsumfun(countsMat)
  if (countsFromAbundance == "lengthScaledTPM") {
    newCounts <- abundanceMat * rowMeans(lengthMat)
  } else if (countsFromAbundance == "scaledTPM") {
    newCounts <- abundanceMat
  } else {
    stop("expecting 'lengthScaledTPM' or 'scaledTPM'")
  }
  newSum <- colsumfun(newCounts)
  if (sparse) {
    countsMat <- Matrix::t(Matrix::t(newCounts) * (countsSum/newSum))
  } else {
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  countsMat
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
          # otherwise use the geometric mean of the lengths from the other samples
          idx <- is.nan(lengthMat[i,])
          lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
        }
    }
  }
  lengthMat
}

cleanTx2Gene <- function(tx2gene) {
  colnames(tx2gene) <- c("tx","gene")
  if (any(duplicated(tx2gene$tx))) {
    message("removing duplicated transcript rows from tx2gene")
    tx2gene <- tx2gene[!duplicated(tx2gene$tx),]
  }
  tx2gene$gene <- factor(tx2gene$gene)
  tx2gene$tx <- factor(tx2gene$tx)
  tx2gene
}

medianLengthOverIsoform <- function(length, tx2gene, ignoreTxVersion, ignoreAfterBar) {
  txId <- rownames(length)
  if (ignoreTxVersion) {
    txId <- sub("\\..*", "", txId)
  } else if (ignoreAfterBar) {
    txId <- sub("\\|.*", "", txId)
  }
  tx2gene <- cleanTx2Gene(tx2gene)
  stopifnot(all(txId %in% tx2gene$tx))
  tx2gene <- tx2gene[match(txId, tx2gene$tx),]
  # average the lengths
  ave.len <- rowMeans(length)
  # median over isoforms
  med.len <- tapply(ave.len, tx2gene$gene, median)
  one.sample <- med.len[match(tx2gene$gene, names(med.len))]
  matrix(rep(one.sample, ncol(length)),
         ncol=ncol(length), dimnames=dimnames(length))
}

summarizeToGeneRaw <- function(txmeta_obj, tx2gene){
  object <- assays(txmeta_obj)
  ignoreTxVersion <- FALSE

  # unpack matrices from list for cleaner code
  abundanceMatTx <- object$abundance
  countsMatTx <- object$counts
  lengthMatTx <- object$length

  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))

  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user

  tx2gene <- cleanTx2Gene(tx2gene)    
    
  # if none of the rownames of the matrices (txId) are
  # in the tx2gene table something is wrong
  if (!any(txId %in% tx2gene$tx)) {
    txFromFile <- paste0("Example IDs (file): [", paste(head(txId,3),collapse=", "),", ...]")
    txFromTable <- paste0("Example IDs (tx2gene): [", paste(head(tx2gene$tx,3),collapse=", "),", ...]")
    stop(paste0("
  None of the transcripts in the quantification files are present
  in the first column of tx2gene."))
  }

  # remove transcripts (and genes) not in the rownames of matrices
  tx2gene <- tx2gene[tx2gene$tx %in% txId,]
  tx2gene$gene <- droplevels(tx2gene$gene)
  ntxmissing <- sum(!txId %in% tx2gene$tx)
  if (ntxmissing > 0) message("transcripts missing from tx2gene: ", ntxmissing)

  # subset to transcripts in the tx2gene table
  sub.idx <- txId %in% tx2gene$tx
  abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
  countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
  lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
  txId <- txId[sub.idx]

  # now create a vector of geneId which aligns to the matrices
  geneId <- tx2gene$gene[match(txId, tx2gene$tx)]


  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- rowsum(abundanceMatTx, geneId)
  message("summarizing counts")
  countsMat <- rowsum(countsMatTx, geneId)
  message("summarizing length")

  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- rowsum(abundanceMatTx * lengthMatTx, geneId)
  lengthMat <- weightedLength / abundanceMat   

  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)

  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))

  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)

  out <- list(counts=countsMat,
                  abundance=abundanceMat,
                  length=lengthMat)
  return(out)
}

########################################################################
# load tximeta object 
########################################################################

if (assembly == "gencodev45") {
  load("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_gencodev45.RData")
  gtf="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"
} else if (assembly == "veigaNorm") {
  load("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_Veiga_ESPRESSO_normal_filtered.RData")
  gtf="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_normal_filtered/ESPRESSO_normal_filtered.gtf"
} else if (assembly == "comb") {
  load("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_combined_veigaNorm_gencodev45.RData")
  gtf="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_normal_filtered.gtf"
} else {
  stop("Incorrect assembly provided. Terminating...")
}

# merge with sample attribute information
sample_dat <- data.frame(full_gtex_id=colnames(se))
sample_dat$gtex_id <- sub("^([^\\-]+\\-[^\\-]+).*", "\\1", sample_dat$full_gtex_id)
sample_attrib <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
sample_dat_new <- merge(x=sample_dat,y=sample_attrib,by.x="full_gtex_id",by.y="SAMPID")
rownames(sample_dat_new) <- sample_dat_new$full_gtex_id
sample_dat_new <- sample_dat_new[sample_dat$full_gtex_id,]
sum(rownames(sample_dat_new)==colnames(se)) # check, 440, good
sample_dat <- sample_dat_new

total_counts <- colSums(assays(se)$counts)
ntx <- nrow(se)
ntx_nonzero <- sum(rowSums(assays(se)$counts)==0)
save(total_counts,ntx, ntx_nonzero,file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/tximeta/GTEx_",assembly,"_persample_read_counts.Rdata"))

# remove low quality samples not intended for QTL
rm <- which(sample_dat$SMAFRZE=="EXCLUDE")
se <- se[,-rm]
sample_dat <- sample_dat[-rm,]
length(unique(sample_dat$gtex_id))==nrow(sample_dat) # true, n=424
sum(sample_dat$full_gtex_id==colnames(se)) # 424, good

# remove transcripts on chr Y & M & X
chr <- as.character(seqnames(rowRanges(se)))
rm <- which(chr %in% c("chrY", "chrM","chrX"))
se <- se[-rm, ]

od <- cs$annotation$Overdispersion
names(od) = rownames(cs$counts)
save(od,file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/tximeta/GTEx_",assembly,"_OD.Rdata"))
od <- od[-rm]

print("dimention of se object:")
dim(se)

# change id to short gtex id
colnames(se) <- sample_dat$gtex_id # there are no repeated short gtex ids

########################################################################
# length scale TPM values if asked
########################################################################
idx_tpm = which(names(assays(se))=="abundance") # abundance is TPM
idx_counts = which(names(assays(se))=="counts")
idx_length = which(names(assays(se))=="length")

if(lstpm==1){
  print("Calculating counts as length-scaled TPM...")
  assays(se)[[idx_counts]] <- makeCountsFromAbundance(countsMat=assays(se)[[idx_counts]],
   abundanceMat=assays(se)[[idx_tpm]],
   lengthMat=assays(se)[[idx_length]],
   countsFromAbundance="lengthScaledTPM")
  }

# counts assay now lengthScaledTPM

########################################################################
# remove lowly expressed transcripts
########################################################################
www = which(rowMeans(assays(se)[[idx_tpm]] > 0.1) > 0.25)
tx_raw = se[www,]
od = od[www]

print("Number of tx meeting TPM threshold:")
dim(tx_raw)

########################################################################
# aggregate tx to gene level & remove lowly expressed genes
########################################################################
gtf <- import(gtf)
df <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])
tx2gene <- distinct(df)

gene_raw <- summarizeToGeneRaw(tx_raw,tx2gene)

www = which(rowMeans(gene_raw[[idx_tpm]] > 0.1) > 0.25)
gene_raw[[1]] = gene_raw[[1]][www,]
gene_raw[[2]] = gene_raw[[2]][www,]
gene_raw[[3]] = gene_raw[[3]][www,]

print("Number of genes meeting TPM threshold:")
print(dim(gene_raw[[1]]))

########################################################################
# ODS scaling
########################################################################

# scale transcript-level counts by overdispersion
tx_counts <- assays(tx_raw)[[idx_counts]]

if(do_ods==1){
  print("Performing ODS scaling...")
  tx_scaled <- sweep(tx_counts, 1, od, "/")
  }else{
    tx_scaled <- tx_counts
  }

########################################################################
# TMM normalization
########################################################################

if(do_tmm==1){
  print("Performing TMM...")
  dge <- DGEList(counts = tx_scaled)
  dge <- calcNormFactors(dge, method = "TMM")
  tx_tmm <- cpm(dge, normalized.lib.sizes = TRUE) # normalized counts in CPM
  #sum(colnames(tx_tmm)==sample_dat$barcode_short)
  #colnames(tx_tmm) <- sample_dat$barcode

  dge <- DGEList(counts = gene_raw[[idx_counts]])
  dge <- calcNormFactors(dge, method = "TMM")
  gene_tmm <- cpm(dge, normalized.lib.sizes = TRUE) # normalized counts in CPM
  #sum(colnames(gene_tmm)==sample_dat$barcode_short)
  #colnames(gene_tmm) <- sample_dat$barcode

  }else{
    tx_tmm <- tx_scaled
    gene_tmm <- gene_raw[[idx_counts]]
  }

########################################################################
# do VST
########################################################################

if(do_vst==1){
  print("Performing VST...")
  tx_cts = DESeq2::vst(round(tx_tmm))
  gene_cts = DESeq2::vst(round(gene_tmm))
  }else{
    tx_cts = tx_tmm
    gene_cts = gene_tmm
  }

########################################################################
# outlier removal
########################################################################

require(WGCNA)
outliers_gene = c()

normadj <- adjacency(gene_cts,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity   #Extract connectivity of each sample
Z.C <- (C-mean(C))/sqrt(var(C))   #Convert to Z-score
outliers <- (Z.C < -3)
print(names(which(outliers == T)))
outliers_gene = c(outliers_gene,
                  names(which(outliers == T)))

outliers <- unique(outliers_gene)

print("Number of outliers detected at gene-level:")
length(outliers)

if(length(outliers)>0){
  rm <- which(sample_dat$gtex_id %in% outliers)
  
  # remove outliers
  tx_cts = tx_cts[,-rm]
  gene_cts = gene_cts[,-rm]
  sample_dat <- sample_dat[-rm,]
  tx_raw = tx_raw[,-rm]
  gene_raw[[1]] = gene_raw[[1]][,-rm]
  gene_raw[[2]] = gene_raw[[1]][,-rm]
  gene_raw[[3]] = gene_raw[[1]][,-rm]

}

print("Samples remaining after outlier removal:")
table(nrow(sample_dat))

########################################################################
# batch correction with ComBat
########################################################################

# THIS IS NOT CORRECT FOR GTEx
# WILL NEED TO FIX
# if(do_bc==1){
#   print("Performing batch correction with ComBat:")
#   sample_dat$TSS_grouped <- ifelse(table(sample_dat$TSS)[sample_dat$TSS] < 3, 
#                                  "Other", 
#                                  sample_dat$TSS)
#   tx_cts <- ComBat(dat = tx_cts,
#     batch = sample_dat$TSS_grouped)

#   gene_cts <- ComBat(dat = gene_cts,
#     batch = sample_dat$TSS_grouped)
# }

########################################################################
# PCA
########################################################################

## PCA on adjusted transcript counts
data_for_pca <- t(tx_cts)
pca_tx <- prcomp(data_for_pca, center = TRUE, scale. = F)

## PCA on adjusted gene counts
data_for_pca <- t(gene_cts)
pca_gene <- prcomp(data_for_pca, center = TRUE, scale. = F)

########################################################################
# write out files
########################################################################

prefix <- paste0(assembly,"_ODS",do_ods,"_TMM",do_tmm,"_VST",do_vst,"_lstpm",lstpm,"_doBC",do_bc)

# save PCA data
save(pca_tx, pca_gene, file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/pca/GTEx_expr_pca_dat_",prefix,".RData"))

# save tx-level data
tx_bed = as.data.frame(rowRanges(tx_raw))
tx_bed$gene.id <- unlist(tx_bed$gene_id)
tx_bed <- tx_bed[,c("seqnames","start","end","tx_name","gene.id","strand")]
tx_bed = cbind(tx_bed,tx_cts)

colnames(tx_bed)[1:5] = c('#Chr','start','end','pid','gid')

print("Ending dimension of tx bed file:")
dim(tx_bed)

tx_bed = tx_bed[order(tx_bed$`#Chr`,
                          tx_bed$start),]
tx_bed$strand <- "+"

require(data.table)
fwrite(tx_bed,paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_tx_exp_',prefix,'.bed'),
       sep='\t',col.names=T,row.names=F)

# save gene-level data

# calculate gene ranges
tx_bed = as.data.frame(rowRanges(tx_raw))
tx_bed$gene.id <- unlist(tx_bed$gene_id)
tx_bed_summary <- tx_bed %>%
  group_by(gene.id) %>%
  summarise(
    min_start = min(start),
    max_end = max(end),
    chr = unique(seqnames)[1]
  )

gene_bed <- data.frame(gene=row.names(gene_raw[[1]]))
gene_bed <- merge(x=gene_bed,y=tx_bed_summary,by.x="gene",by.y="gene.id",all.x=T,all.y=F)
row.names(gene_bed) <- gene_bed$gene
gene_bed <- gene_bed[row.names(gene_cts),]
sum(row.names(gene_bed)==row.names(gene_cts))==nrow(gene_cts) # good

gene_bed <- gene_bed[,c("chr","min_start","max_end","gene","gene")]
gene_bed$strand <- "+"
gene_bed = cbind(gene_bed,gene_cts)

colnames(gene_bed)[1:5] = c('#Chr','start','end','pid','gid')

print("Ending dimension of tx bed file:")
dim(gene_bed)

gene_bed = gene_bed[order(gene_bed$`#Chr`,
                          gene_bed$start),]

fwrite(gene_bed,paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_gene_exp_',prefix,'.bed'),
       sep='\t',col.names=T,row.names=F)
