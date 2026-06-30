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
library(isotwas)
library(data.table)
library(bigsnpr)
library(vroom)
library(glmnet)
library(rrBLUP)
library(susieR)
library(limma)
library(rtracklayer)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly = TRUE)
exp_file <- as.character(args[1])
gen_file <- as.character(args[2])
cov_file <- as.character(args[3])
label <- as.character(args[4])
gtf_file <- as.character(args[5])
chr <- as.numeric(args[6])
index <- as.numeric(args[7])

# exp_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/TCGA_tx_exp_comb_ODS1_TMM1_VST1_lstpm0_doBC0_tumor.bed.gz"
# gen_file="/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03/autosome_unrelated_full_ids"
# cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/TCGA_cov_qtltools_tx_tumor_comb_ODS1_TMM1_VST1_lstpm0_doBC0_k15_HC_genpc5.txt"
# label="comb_k15"
# gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_tumor_filtered.gtf"
# chr=1
# index=1

########################################################################
# begin code
########################################################################

annot <- strsplit(label,"_")[[1]][1]

###ISOFORM EXPRESSION
exp = vroom(exp_file)
exp = as.data.frame(exp)

###COVARIATE
cov <- data.frame(fread(cov_file,header=F))
colnames(cov) <- cov[1,]
rownames(cov) <- cov[,1]
#cov <- cov[-1,-1]
cov <- cov[-1,]
shared_samples <- intersect(colnames(cov)[2:ncol(cov)], colnames(exp)[7:ncol(exp)])

#Subset each dataframe to only shared samples, and order columns identically
cov_sub <- cov[, c(colnames(cov)[1], shared_samples)]  # keep 1st column (e.g., ID/feature row)
exp_sub <- exp[, c(colnames(exp)[1:6], shared_samples)] # keep first 6 columns of BED

### COVARIATE RESIDUALIZE
cov_df = data.frame(t(cov_sub))
colnames(cov_df) = cov_sub$id
cov_df = cov_df[-1,]
cov_df = as.data.frame(apply(cov_df,2,as.numeric))
#cov_df$sex = as.character(cov_sub[cov_sub$id == 'sex',-1])
mm = model.matrix(~.,data = cov_df)
rm(cov,cov_sub,cov_df,exp)
exp_sub[,7:ncol(exp_sub)] = removeBatchEffect(as.matrix(exp_sub[,7:ncol(exp_sub)]),
                                              covariates = mm[,-1])


###GENE ID

# if novel transcripts mapped to multiple genes, create multiple rows for each one
if(annot =="veigaTum"){

  exp_sub <- exp_sub[-which(is.na(exp_sub$gid)),] # remove novel transcripts not assigned to a gene

  multiple_genes <- grep(",",exp_sub$gid)
  exp_fine <- exp_sub[-multiple_genes,]
  exp_mult <- exp_sub[multiple_genes,]

  exp_extended <- {}
  for(i in 1:nrow(exp_mult)){
    tmp <- exp_mult[i,]
    gene_ids <- strsplit(tmp$gid,",")[[1]]

    new <- tmp[rep(1,length(gene_ids)),]
    new$gid <- gene_ids
    exp_extended <- rbind(exp_extended,new)
  }

  exp_sub <- rbind(exp_fine,exp_extended)
  exp_sub <- exp_sub[order(exp_sub$start,decreasing=F),]
  exp_sub <- exp_sub[order(exp_sub$`#Chr`),]
}

bed_gene = unique(exp_sub$gid)

gtf <- import(gtf_file)
gtf <- data.frame(gtf)
#gtf <- gtf[gtf$type == "gene",]
gtf <- gtf[gtf$gene_id %in% bed_gene,]
gtf <- gtf[gtf$seqnames==paste0("chr",chr),]

gene_ids <- unique(gtf$gene_id)
gene_ids <- sort(gene_ids)
n_bins <- 100
gene_bins <- split(gene_ids, cut(seq_along(gene_ids), breaks = n_bins, labels = FALSE))
gene_ids <- gene_bins[[index]]
gtf_batch <- gtf[gtf$gene_id %in% gene_ids,]
exp_sub = subset(exp_sub,gid %in% gene_ids)

### READ (PARTIAL) GENOTYPE DATA
pvar_path <- paste0(gen_file, ".pvar")
pgen_path <- paste0(gen_file, ".pgen")
psam_path <- paste0(gen_file, ".psam")
pvar <- pgenlibr::NewPvar(pvar_path) 
pgen <- pgenlibr::NewPgen(pgen_path, pvar=pvar) # Check the number of variants and samples. pgenlibr::GetVariantCt(pgen) pgenlibr::GetRawSampleCt(pgen)
pvar_mat <- fread(paste0(gen_file,".pvar"))
pvar_idx <- which(pvar_mat$`#CHROM`==chr)

pgenlibr::GetVariantCt(pgen) 
pgenlibr::GetRawSampleCt(pgen)

geno_mat <- pgenlibr::ReadList(pgen,pvar_idx)
snp_info <- pvar_mat[pvar_idx,]

colnames(geno_mat) <- snp_info$ID

ss <- fread("/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_GRCh38_overallBC_ambig_strand_removed.tsv.gz")
ss <- ss[ss$CHR ==chr,]
ss$ID <- paste(ss$CHR,ss$BP,ss$A1,ss$A2,sep=":")
ss$ID2 <- paste(ss$CHR,ss$BP,ss$A2,ss$A1,sep=":")

keep <- which(snp_info$ID %in% ss$ID)
keep2 <- which(snp_info$ID %in% ss$ID2)
keep <- unique(c(keep,keep2))
snp_info <- snp_info[keep,]
geno_mat <- geno_mat[,keep]

sample_ids <- data.frame(fread(psam_path))
sample_ids <- substr(sample_ids[,1],1,12)
rownames(geno_mat) <- sample_ids

#snps$genotypes <- snp_fastImputeSimple(snps$genotypes, method = "mean2")

runTCGA_isotwas = function(gene){
  #gene = "ENSG00000000419.14"
  exp_current = exp_sub[exp_sub$gid %in% gene,]
  
  if (nrow(exp_current) > 0){
    
    exp_mat = as.matrix(t(as.matrix(exp_current[,7:ncol(exp_current)])))
    colnames(exp_mat) = exp_current$pid
    #rownames(exp_mat) = colnames(exp_current[7:ncol(exp_current)])
    chr = unique(exp_current$`#Chr`)[1]
    start = max(1,min(exp_current$start) - 1e6)
    end = max(exp_current$start) + 1e6

    snp_list = snp_info[snp_info$POS <= end & snp_info$POS >= start,]
    snp_list = snp_list$ID
    # snp_list = snps$map$marker.ID[snps$map$chromosome == chr &
    #                                 snps$map$physical.pos < end &
    #                                 snps$map$physical.pos > start]
    
    if (length(snp_list) > 1){

      snp_current = geno_mat[,snp_list]
      # snp_current = snp_attach(subset(snps,
      #                                 ind.col = which(snps$map$marker.ID %in%
      #                                                   snp_list)))
      
      ### Train isoTWAS model 
  
      snpMat = as.matrix(snp_current)
      snpMat = snpMat[rownames(exp_mat),]
      
      if (ncol(exp_mat) > 1){
        
        m.isot = compute_isotwas(X = snpMat,
                                 Y = exp_mat,
                                 Y.rep = exp_mat,
                                 R = 1,
                                 gene_exp = NULL,
                                 id = rownames(snpMat),
                                 omega_est = 'replicates',
                                 omega_nlambda = 5,
                                 method = c('mrce_lasso',
                                            'multi_enet',
                                            'univariate',
                                            'mvsusie'),
                                 predict_nlambda = 5,
                                 family = 'gaussian',
                                 scale = F,
                                 alpha = 0.5,	
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = colnames(exp_mat),
                                 seed = 1218,
                                 run_all = F,
                                 return_all = T,
                                 tol.in = 1e-3,
                                 coverage = .9)
        tx2gene = m.isot$tx2gene_coef
        m.isot = m.isot$isotwas_mod
        
        model.df = data.frame(Feature = NA,
                              SNP = NA,
                              Chromosome = NA,
                              Position = NA,
                              Build = NA,
                              ALT = NA,
                              REF = NA,
                              Weight = NA,
                              R2 = NA)
        for (j in 1:length(m.isot$Model)){
          
          aaa = m.isot$Model[[j]]
          if (class(aaa) == 'data.frame'){
            this_model = aaa
            
            if (nrow(this_model) > 0){
              this_model$Feature = exp_current$pid[1]
              model.df = rbind(model.df,this_model)
            }
            
          } else {
            if (nrow(aaa$Model) == 0) next 
            colnames(aaa$Model) = c('marker.ID','Weight')
            aaa$Model = merge(aaa$Model,
                              snp_info,
                              by.x = 'marker.ID',
                              by.y="ID")
            aaa$Model = aaa$Model[,c('marker.ID','Weight',
                                     '#CHROM','POS',
                                     'REF','ALT')]
            colnames(aaa$Model) = c('SNP',
                                    'Weight',
                                    'Chromosome',
                                    'Position',
                                    'REF',
                                    'ALT')
            aaa$Model$R2 = aaa$R2
            aaa$Model$Build = 'hg38'
            aaa$Model$Feature = colnames(exp_mat)[j]
            aaa$Model = aaa$Model[,c('Feature','SNP','Chromosome','Position','Build',
                                     'ALT','REF','Weight','R2')]
            model.df = rbind(model.df,aaa$Model)
            
            
          }
        }
        
      } else {
        
        exp_mat = as.matrix(t(as.matrix(exp_current[,7:ncol(exp_current)])))
        colnames(exp_mat) = exp_current$pid
        #rownames(exp_mat) = colnames(exp_current[7:ncol(exp_current)])
        chr = unique(exp_current$`#Chr`)[1]
        start = max(1,min(exp_current$start) - 1e6)
        end = max(exp_current$start) + 1e6
        snp_list = snp_info[snp_info$POS < end & snp_info$POS > start,]
        snp_list = snp_list$ID
        
        snp_current = geno_mat[,snp_list]
        snpMat = as.matrix(snp_current)
        snpMat = snpMat[rownames(exp_mat),]
        
        enet = univariate_elasticnet(X = snpMat,
                                     Y = exp_mat,
                                     Omega = NA,
                                     family = 'gaussian',
                                     scale = F,
                                     alpha = 0.5,
                                     nfolds = 5,
                                     verbose = F,
                                     par = F,
                                     n.cores = NULL,
                                     tx_names = NULL,
                                     seed = NULL)
        
        blup = univariate_blup(X = snpMat,
                               Y = exp_mat,
                               Omega = NA,
                               scale = F,
                               #alpha = 0.5,
                               nfolds = 5,
                               verbose = F,
                               par = F,
                               n.cores = NULL,
                               tx_names = NULL,
                               seed = NULL)
        
        susie = univariate_susie(X = snpMat,
                                 Y = exp_mat,
                                 Omega = NA,
                                 scale = F,
                                 alpha = 0.5,
                                 nfolds = 5,
                                 verbose = F,
                                 par = F,
                                 n.cores = NULL,
                                 tx_names = NULL,
                                 seed = NULL)
        
        r2 = unlist(c(enet[[1]]$R2,
                      blup[[1]]$R2,
                      susie[[1]]$R2))
        r2.min = which(r2 == max(r2))
        if (r2.min == 1){
          last.mod = enet
          if (nrow(enet[[1]]$Model) == 0){
            last.mod = susie
            r2.min = 3
          }
        }
        
        if (r2.min == 2){
          last.mod = blup
        }
        
        if (r2.min == 3){
          last.mod = susie
        }
        
        gene.model = last.mod[[1]]$Model
        gene.model$Feature = last.mod[[1]]$Transcript
        gene.model$R2 = last.mod[[1]]$R2
        gene.model$marker.ID = gene.model$SNP
        gene.model = merge(gene.model,snp_info,by.x='marker.ID',by.y="ID")
        colnames(gene.model) = c('marker.ID','SNP',
                                 'Weight','Feature','R2','Chromosome',
                                 'Position',
                                 'REF','ALT')
        gene.model$Build = 'hg38'
        gene.model = gene.model[,c('Feature','SNP','Chromosome','Position','Build',
                                   'ALT','REF','Weight','R2')]
        model.df = gene.model
        
      }
      
      model.df = model.df[!is.na(model.df$R2),]
      #ld.cor = Matrix::Matrix(snp_cor(snp_current))
      ld.cor <- cor(snp_current, use = "pairwise.complete.obs")
      ld.cor <- Matrix::Matrix(ld.cor)
      #colnames(ld.cor) = rownames(ld.cor) = snp_current$map$marker.ID
      
      out_folder = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged')
      out_folder_isoTWAS = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/',label,'_isoTWAS')
      out_folder_ld = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/',label,'_LDMatrix')
      out_folder_doneFile = '/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/doneFile'
      dir.create(out_folder,recursive = T,showWarnings = F)
      dir.create(out_folder_isoTWAS,recursive = T,showWarnings = F)
      dir.create(out_folder_ld,recursive = T,showWarnings = F)
      dir.create(out_folder_doneFile,recursive = T,showWarnings = F)
      
      saveRDS(model.df,file.path(out_folder_isoTWAS,
                                 paste0(gene,'_isoTWAS.RDS')))
      saveRDS(ld.cor,
              file.path(out_folder_ld,
                        paste0(gene,'_LDMatrix.RDS')))
    }
  }
}


doneFile = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/doneFile/done_models_', label, '.tsv')

if (file.exists(doneFile)) {
  done = fread(doneFile)
} else {
  done = data.table(Gene = character(), Done = integer())
}

for (gene in gene_ids) {
  print(gene)
  if (!(gene %in% done$Gene)) {
    tryCatch({
      runTCGA_isotwas(gene)
    }, error = function(e) {
      message("Gene did not fit: ", gene)
    })

    done_df = data.frame(Gene = gene, Done = 1)
    fwrite(done_df, doneFile, append = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
}