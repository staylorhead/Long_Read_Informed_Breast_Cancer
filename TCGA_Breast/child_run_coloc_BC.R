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
library(coloc)
library(data.table)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
gwas_file <- as.character(args[1])
eqtl_file <- as.character(args[2])
pheno <- as.character(args[3])
chunk <- as.integer(args[4])
level <- as.character(args[5])
annot <- as.character(args[6])

####################################################################################
# begin
####################################################################################
setwd("/rsrch5/home/epi/sthead/isoqtl_lr_breast")

eqtl <- fread(eqtl_file)
header <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/nominal_header_0.txt")  # update as needed
setnames(eqtl, names(header))

gwas <- fread(gwas_file)
gwas[, var_id := paste(CHR, BP, A1, A2, sep = "_")]

# read in psam file to change rsid to var id
pvar <- fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03/allchr.pvar")
#pvar[, chr := paste0("chr", `#CHROM`)]
pvar[, var_id_new := paste(`#CHROM`, POS, REF, ALT, sep = "_")]
eqtl_annot <- merge(eqtl, pvar[, .(ID, var_id_new)], by.x = "var_id", by.y = "ID", all.x = TRUE)
eqtl <- eqtl_annot

# Clean and harmonize
eqtl <- eqtl[!grepl("^HGSV", var_id)]
eqtl[, c("ref_allele", "alt_allele") := tstrsplit(var_id_new, "_", keep = 3:4)]
eqtl[, var_chr_nochr := sub("^chr", "", var_chr)]
eqtl[, flipped_var_id := paste(var_chr_nochr, var_from, alt_allele, ref_allele, sep = "_")]
eqtl[, var_id_in_gwas := var_id_new %in% gwas$var_id]
eqtl[, flipped_var_id_in_gwas := flipped_var_id %in% gwas$var_id]
eqtl[, slope_var := slope_se^2]

# Flip alleles where needed
eqtl[var_id_in_gwas == FALSE & flipped_var_id_in_gwas == TRUE,
     `:=`(var_id_new = flipped_var_id,
          ref_allele = alt_allele,
          alt_allele = ref_allele,
          slope = -slope)]

# Get phe_ids to iterate through
phe_ids <- unique(eqtl$phe_id)

# Initialize results list
coloc_results <- list()

gwas[, Z := BETA / SE]
gwas[, MAF := pmin(FRQ, 1 - FRQ)]

if (!"VAR" %in% names(gwas)) {
  gwas[, VAR := SE^2]
}

if(pheno=='overallBC'){
  s=133384/ (133384+113789)
  N=(133384+113789)
  }else if(pheno=="lumA"){
  s=106278/(106278+91477)
  N=(106278+91477)
    }else if(pheno=="lumB"){
      s=106278/(106278+91477)
      N=(106278+91477)
    }else if (pheno=="lumBHER2neg"){
      s=106278/(106278+91477)
      N=(106278+91477)
      }else if(pheno=="HER2enriched"){
        s=106278/(106278+91477)
        N=(106278+91477)
      }else if(pheno=="TN"){
        s = 16592 / (16592 + 103901)
        N=(16592 + 103901)
}

# Iterate over phe_id
for (phe in phe_ids) {

  eqtl_phe <- eqtl[phe_id == phe]
  eqtl_phe <- eqtl_phe[!eqtl_phe$slope_se==0,]
  if (nrow(eqtl_phe)==0) next

  phe_start <- max(0, eqtl_phe$phe_from[1] - 1e6)
  phe_end   <- eqtl_phe$phe_to[1] + 1e6
  phe_chr   <- as.integer(eqtl_phe$var_chr_nochr[1])

  # Subset GWAS region
  gwas_phe <- gwas[CHR == phe_chr & BP >= phe_start & BP <= phe_end]
  if (nrow(gwas_phe)==0) next

  # Compute MAF and var(beta)
  gwas_phe[, var_g := 2 * FRQ * (1 - FRQ)]
  #gwas_phe[, beta := Z / sqrt(N * var_g)]
  #gwas_phe[, var_beta := 1 / (N * var_g)]

  # Ensure SNPs match
  overlap_snps <- intersect(eqtl_phe$var_id_new, gwas_phe$var_id)
  if (length(overlap_snps) < 2) next

  # Prepare coloc input
  eqtlList <- list(
    snp = eqtl_phe$var_id_new,
    beta = eqtl_phe$slope,
    varbeta = eqtl_phe$slope_var,
    N = 1083,  # or adjust dynamically
    type = "quant",
    pvalues=eqtl_phe$nom_pval,
    sdY = 1
  )

  # gwasList <- list(
  #   snp = gwas_phe$var_id,
  #   MAF = gwas_phe$MAF,
  #   Z = gwas_phe$Z,
  #   N = N,
  #   s=s,
  #   pvalues = gwas_phe$P,
  #   type = "cc"
  # )

    gwasList <- list(
    snp = gwas_phe$var_id,
    MAF = gwas_phe$MAF,
    beta = gwas_phe$BETA,
    varbeta = gwas_phe$VAR,
    s=s,
    type = "cc"
  )

  res <- coloc.abf(dataset1 = eqtlList, dataset2 = gwasList)
  coloc_results[[phe]] <- res
}

dir_path <- "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc/"

if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

# Save output (e.g., per chunk)
save(coloc_results, file = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc/TCGA_",annot,"_",pheno,"_",level,"_results_chunk", chunk, ".RData"))







