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
library(data.table)
library(dplyr)
library(pgenlibr)
library(susieR)
library(Matrix)

####################################################################################
# parse arguments
####################################################################################
args <- commandArgs(trailingOnly = TRUE)
level <- as.character(args[1]) # gene or tx
annot <- as.character(args[2])
k <- as.numeric(args[3])
nperm <- as.numeric(args[4])
maf_filt <- as.numeric(args[5])
rm_cov <- as.numeric(args[6])
nsets <- as.numeric(args[7])
i <- as.numeric(args[8])
first_numeric_cov <- as.numeric(args[9])
L <- as.numeric(args[10])

# level="tx"
# annot="veigaNorm"
# k=15
# nperm=1000
# maf_filt=0.01
# rm_cov=1
# nsets=100
# i=1
# first_numeric_cov=2
# L=5

####################################################################################
# begin
####################################################################################

# Univariate SuSiE Analysis for a Single Transcript
# =======================================================
# This script:
# 1. Extracts gene data using PLINK to subset genotype data
# 2. Selects a single transcript for analysis
# 3. Runs univariate SuSiE to identify SNPs affecting the transcript
# 4. Loops througha all transripts or genes in a set

# File paths
qtlSummaryStats = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools/GTEx_fibro_cond_",annot,"_ODS1_TMM1_VST1_lstpm0_doBC0_k",k,"_HC_genpc5/",level,"_qtls_cond_nperm",nperm,"_normalized_full.txt.gz")
genFile = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/geno/GTEx_MAF0.01_passQC/autosome_unrelated_sorted_chr_labels.vcf.gz")
expressionBed = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_fibroblast_",level,"_exp_",annot,"_ODS1_TMM1_VST1_lstpm0_doBC0.bed.gz")
covFile = paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_fibro_cov_qtltools_",level,"_",annot,"_ODS1_TMM1_VST1_lstpm0_doBC0_k",k,"_HC_genpc5.txt")
condHeader = fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/conditional_header_0.txt")

# Extract column names from covFile for sample IDs
covData = fread(covFile)
covNames <- c("sex","age",paste0("hcp",1:k),paste0("PC",1:5))

rownames(covData) <- covNames
covDataGeno <- covData[c(grep("sex",rownames(covData)),grep("PC",rownames(covData))),]

masterListIDs = names(covData)[-1]  # all column names except the first

# Read data
qtlSS = fread(qtlSummaryStats)
colnames(qtlSS) <- colnames(condHeader)

sigTx = unique(subset(qtlSS, bwd_sig == 1)$phe_id)

# subset to set of 100
chunks <- split(sigTx, cut(seq_along(sigTx), breaks = nsets, labels = FALSE))
# get the i-th chunk for this job
sigTx <- chunks[[i]]

tx2Gene = fread(expressionBed, select = 1:6)
colnames(tx2Gene)[1:6] <- c("chr", "start", "end", "pid", "gid", "strand")
tx2Gene = subset(tx2Gene, pid %in% sigTx)

#################################################
# Function to remove covariate effects
#################################################

# this function regresses out covariates froms genotype 
# remove.covariate.effects <- function (X, Z) {
#   # include the intercept term
#   if (any(Z[,1]!=1)) Z = cbind(1, Z)
#   A   <- forceSymmetric(crossprod(Z))
#   # SZy <- as.vector(solve(A,c(y %*% Z)))
#   SZX <- as.matrix(solve(A,t(Z) %*% X))
#   # y <- y - c(Z %*% SZy)
#   X <- X - Z %*% SZX
#   return(list(X = X,SZX = SZX))
# }
remove.covariate.effects <- function(X, Z) {
  # Ensure numeric matrix
  Z <- as.matrix(data.frame(lapply(Z, function(x) as.numeric(as.character(x)))))
  # Drop any all-NA or constant columns
  keep <- apply(Z, 2, function(v) sum(!is.na(v)) > 0 && length(unique(v)) > 1)
  Z <- Z[, keep, drop = FALSE]
  # Add intercept
  Z <- cbind(1, Z)
  
  # Residualize
  A <- crossprod(Z)
  SZX <- solve(A, crossprod(Z, X))
  X <- X - Z %*% SZX
  return(list(X = X, SZX = SZX))
}

#################################################
# Function to extract data for a specific transcript
#################################################

extract_transcript_data <- function(transcript_id) {
  # Check if transcript exists
  if (!(transcript_id %in% tx2Gene$pid)) {
    stop(paste("Transcript", transcript_id, "not found in dataset."))
  }
  
  # Get transcript information
  transcript_info <- subset(tx2Gene, pid == transcript_id)
  target_gene <- transcript_info$gid[1]
  
  cat("Selected transcript:", transcript_id, "from gene:", target_gene, "\n")
  
  # Get gene information for region extraction
  min_start <- transcript_info$start[1]
  gene_chr <- transcript_info$chr[1]
  
  # Strip 'chr' prefix if present in gene_chr
  gene_chr <- gsub("^chr", "", gene_chr)
  
  # Define the region to extract SNPs (1Mb upstream)
  region_start <- max(0, min_start - 1e6)
  region_end <- min_start + 1e6  # 1Mb downstream as well for broader context
  
  cat("Extracting SNPs from chromosome", gene_chr, "position", region_start, "to", region_end, "\n")
  
  # Create a region string for extracting SNPs
  region_string <- paste0(gene_chr, ":", region_start, "-", region_end)
  
  # Create a temporary prefix for the output files
  temp_prefix <- tempfile()
  #temp_prefix <- tempfile(tmpdir = "/scratch/sthead/tmp")

  # Use PLINK2 to extract the region
  if(!is.na(maf_filt)){

    plink_cmd <- paste0(
    "plink2 --vcf ", genFile, " ",
    "--chr ", gene_chr, " ",
    "--from-bp ", region_start, " ",
    "--to-bp ", region_end, " ",
    "--maf ",maf_filt," ",  # Filter out variants with MAF < threshold
    "--make-pgen ",  # Output in PGEN format
    "--out ", temp_prefix)
    }else{

    plink_cmd <- paste0(
    "plink2 --vcf ", genFile, " ",
    "--chr ", gene_chr, " ",
    "--from-bp ", region_start, " ",
    "--to-bp ", region_end, " ",
    "--make-pgen ",  # Output in PGEN format
    "--out ", temp_prefix)
  }
  
  # Execute PLINK command
  cat("Executing PLINK2 command:\n", plink_cmd, "\n")
  system(plink_cmd, ignore.stdout = T, ignore.stderr = T)
  
  # Check if PLINK command was successful and files were created
  if (!file.exists(paste0(temp_prefix, ".pgen"))) {
    stop("PLINK filtering failed. Check if region contains any SNPs.")
  }
  
  # Read the filtered genotype data using pgenlibr
  pgen <- NewPgen(paste0(temp_prefix, ".pgen"))
  pvar <- paste0(temp_prefix, ".pvar")
  psam <- paste0(temp_prefix, ".psam")
  
  map = fread(pvar)
  fam = fread(psam)
  
  # Read the genotype data
  geno_matrix <- as.matrix(ReadList(pgen, variant_subset=1:nrow(map), meanimpute=T))
  
  # Set proper row and column names
  colnames(geno_matrix) <- map$ID
  rownames(geno_matrix) <- fam$`#IID`
  
  # Restrict rownames to first 12 characters
  rownames(geno_matrix) <- substr(rownames(geno_matrix), 1, 12)
  
  # Keep only rows that match masterListIDs
  geno_matrix <- geno_matrix[rownames(geno_matrix) %in% masterListIDs, ]
  
  # Reorder to match the exact order in masterListIDs
  matched_ids <- intersect(masterListIDs, rownames(geno_matrix))
  geno_matrix <- geno_matrix[match(matched_ids, rownames(geno_matrix)), ]
  
  # Read expression data for the specific transcript
  col_count <- length(fread(expressionBed, nrows = 1))
  col_select <- c(1:6, 7:col_count)
  
  # Read the full expression file
  expression_data <- fread(expressionBed, select = col_select)
  # Rename first 6 columns
  colnames(expression_data)[1:6] <- c("chr", "start", "end", "pid", "gid", "strand")
  
  # Filter for the specific transcript
  transcript_expression <- expression_data[pid == transcript_id]
  
  # Extract just the expression values (all columns except the first 6)
  expression_vector <- as.numeric(as.vector(transcript_expression[, 7:ncol(transcript_expression)]))
  
  # Get the sample IDs
  sample_ids <- colnames(expression_data)[7:ncol(expression_data)]
  
  # Create a named vector of expression values
  names(expression_vector) <- sample_ids
  
  # Filter for samples that match masterListIDs
  expression_vector <- expression_vector[names(expression_vector) %in% masterListIDs]
  
  # Reorder to match the exact order in masterListIDs
  matched_samples <- intersect(masterListIDs, names(expression_vector))
  expression_vector <- expression_vector[match(matched_samples, names(expression_vector))]
  
  # Clean up temporary files and close pgenlibr objects
  cleanup_files <- function() {
    # Close pgenlibr objects
    if (exists("pgen")) ClosePgen(pgen)
    
    # Remove temporary files
    file_extensions <- c(".pgen", ".pvar", ".psam", ".log")
    for (ext in file_extensions) {
      file_path <- paste0(temp_prefix, ext)
      if (file.exists(file_path)) {
        file.remove(file_path)
      }
    }
  }
  
  # Register cleanup to happen even if there's an error
  on.exit(cleanup_files(), add = TRUE)
  
  # Return the results
  return(list(
    transcript_id = transcript_id,
    gene_id = target_gene,
    genotypes = geno_matrix,
    expression = expression_vector
  ))
}

#################################################
# Function to run univariate SuSiE analysis
#################################################

run_univariate_susie <- function(transcript_id, L = 10) {
  # Extract data for the transcript
  cat("Extracting data for transcript:", transcript_id, "\n")
  transcript_data <- extract_transcript_data(transcript_id)
  
  # Get X (genotypes) and y (expression)
  X <- transcript_data$genotypes
  y <- transcript_data$expression

  # Ensure X and y have matching samples
  common_samples <- intersect(intersect(rownames(X), names(y)),masterListIDs)
  X <- X[common_samples, ]
  y <- y[common_samples]

  # Prepare Z (covariates)
  tcov <- data.frame(covData[,-1])
  colnames(tcov) <- masterListIDs
  Z <- t(tcov[, common_samples])
  Z <- as.data.frame(Z)

  # Prepare Z_geno (covariates for genotype data)
  tcov_geno <- data.frame(covDataGeno[,-1])
  colnames(tcov_geno) <- masterListIDs
  Z_geno <- t(tcov_geno[, common_samples])
  Z_geno <- as.data.frame(Z_geno)

  # Convert covariates
  if ((first_numeric_cov - 1) >= 1) {
    Z[, 1:(first_numeric_cov - 1)] <- as.data.frame(lapply(Z[, 1:(first_numeric_cov - 1), drop = FALSE], as.factor))
  }
  Z[, first_numeric_cov:ncol(Z)] <- as.data.frame(lapply(Z[, first_numeric_cov:ncol(Z), drop = FALSE], as.numeric))

  # Convert covariates
  # Convert non-numeric covariates in Z_geno
  if ((first_numeric_cov - 1) >= 1) {
    Z_geno[, 1:(first_numeric_cov - 1)] <- lapply(
      Z_geno[, 1:(first_numeric_cov - 1), drop = FALSE],
      function(x) {
        if (is.factor(x) || is.character(x)) {
          as.numeric(factor(x)) - 1
        } else {
          as.numeric(x)
        }
      }
    )
  }
  Z_geno[, first_numeric_cov:ncol(Z_geno)] <- as.data.frame(lapply(Z_geno[, first_numeric_cov:ncol(Z_geno), drop = FALSE], as.numeric))

  cat("Data dimensions - X:", dim(X), "y:", length(y), "\n")
  
  # Check if there are enough SNPs to proceed
  if (ncol(X) < 10) {
    cat("Too few SNPs in the region. Skipping analysis.\n")
    return(NULL)
  }
  
  # Standardize X
  X_std <- scale(X, center = TRUE, scale = TRUE)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Run univariate SuSiE
  cat("Running univariate SuSiE with L =", L, "effects...\n")
  
  # Removing covariate effects
  if(rm_cov==1){
    out = remove.covariate.effects(X_std, Z_geno)
    X_std <- out$X
    # y <- out$y
    y <- residuals(lm(y ~ ., data = as.data.frame(Z)))
  }
  
  susie_fit <- susie(
    X = X_std,
    y = y,
    L = L,                        # Number of effects to fit
    standardize = FALSE,          # Already standardized
    verbose = TRUE,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE
  )
  
  cat("SuSiE analysis complete.\n")
  
  # Return the results
  return(list(
    transcript_id = transcript_id,
    gene_id = transcript_data$gene_id,
    susie_fit = susie_fit,
    X_snps = colnames(X_std)
  ))
}

#################################################
# Function to summarize SuSiE results
#################################################

summarize_susie_results <- function(result) {
  if (is.null(result)) {
    return("Analysis could not be completed due to insufficient data.")
  }
  
  # Extract SuSiE fit
  susie_fit <- result$susie_fit
  
  # Get the posterior inclusion probabilities (PIPs)
  pips <- susie_fit$pip
  
  # Get credible sets
  cs <- susie_get_cs(susie_fit)
  
  # Create a summary
  summary_text <- paste0(
    "=== Summary of SuSiE results ===\n",
    "Transcript ID: ", result$transcript_id, "\n",
    "Gene ID: ", result$gene_id, "\n",
    "Number of variables with PIP > 0.5: ", sum(pips > 0.5), "\n",
    "Number of credible sets: ", length(cs$cs), "\n\n"
  )
  
  # Add information about top SNPs
  if (length(pips) > 0) {
    summary_text <- paste0(summary_text, "=== Top SNPs by PIP ===\n")
    top_indices <- order(pips, decreasing = TRUE)[1:min(5, length(pips))]
    
    for (i in 1:length(top_indices)) {
      idx <- top_indices[i]
      snp_id <- names(pips)[idx]
      pip_val <- pips[idx]
      summary_text <- paste0(summary_text, i, ". ", snp_id, " (PIP = ", round(pip_val, 4), ")\n")
    }
  }
  
  # Add information about credible sets
  if (length(cs$cs) > 0) {
    summary_text <- paste0(summary_text, "\n=== Credible Sets ===\n")
    
    for (i in 1:length(cs$cs)) {
      cs_snps <- cs$cs[[i]]
      cs_size <- length(cs_snps)
      
      # Get the top SNP in the credible set
      top_snp_idx <- cs_snps[1]
      top_snp_id <- names(pips)[top_snp_idx]
      top_snp_pip <- pips[top_snp_idx]
      
      summary_text <- paste0(
        summary_text, 
        "Credible Set ", i, " (", cs_size, " SNPs)\n",
        "  Top SNP: ", top_snp_id, " (PIP = ", round(top_snp_pip, 4), ")\n"
      )
      
      if (cs_size > 1) {
        summary_text <- paste0(
          summary_text,
          "  Also includes ", cs_size - 1, " additional SNPs\n"
        )
      }
    }
  }
  
  # Add information about estimated effect sizes
  if (!is.null(susie_fit$alpha) && !is.null(susie_fit$mu)) {
    summary_text <- paste0(summary_text, "\n=== Estimated Effect Sizes ===\n")
    
    for (l in 1:length(susie_fit$KL)) {
      # Alpha is the weight for each variable in the lth component
      # Mu is the conditional mean for each variable in the lth component
      weighted_effects <- susie_fit$alpha[, l] * susie_fit$mu[, l]
      
      # Get the SNP with the largest weight in this component
      top_idx <- which.max(abs(weighted_effects))
      top_snp <- names(pips)[top_idx]
      effect_size <- weighted_effects[top_idx]
      
      summary_text <- paste0(
        summary_text,
        "Effect ", l, ": ",
        top_snp, " (Effect Size = ", round(effect_size, 4), ")\n"
      )
    }
  }
  
  return(summary_text)
}

#################################################
# Example usage
#################################################

# target_transcript <- sigTx[1]
# cat("Selected transcript for analysis:", target_transcript, "\n")

# # Define a wrapper for your fine-mapping function
# run_finemap_wrapper <- function(target_transcript) {
#     message("Running fine-mapping for: ", target_transcript)

#     # Call your actual fine-mapping logic here
#     result <- run_univariate_susie(target_transcript, L = 10)
#     return(result)
# }

# # Use lapply to process each feature in this job
# fullres <- lapply(sigTx, run_finemap_wrapper)

fullres <- list()
for (tx in sigTx) {
  cat("Processing", tx, "...\n")
  result <- tryCatch({
    run_univariate_susie(transcript_id=tx,L=L)
  }, error = function(e) {
    cat("Error processing", tx, ":", e$message, "\n")
    return(NULL)  # or return a list with error info
  })
  fullres[[tx]] <- result
}

dir_path <- "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap"

if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

save(fullres,file=paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap/GTEx_fibro_",level,"_",annot,"_k",k,"_nperm",nperm,"_maf_filt",maf_filt,
  "_rmcov",rm_cov,"_",i,"_L",L,".RData"))


