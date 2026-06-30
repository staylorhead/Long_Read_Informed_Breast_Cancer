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
library(isotwas)
library(data.table)
library(bigsnpr)
library(dplyr)
library(Matrix)

####################################################################################
# parse arguments
####################################################################################

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
bin_index <- as.numeric(args[1])
isotwas_folder <- as.character(args[2])
ld_folder <- as.character(args[3])
gwas_file <- as.character(args[4])
label <- as.character(args[5])
trait <- as.character(args[6])

# isotwas_folder="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/gencodev45_k15_isoTWAS"
# ld_folder = "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/gencodev45_k15_LDMatrix"
# gwas_file="/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_GRCh38_overallBC_ambig_strand_removed.tsv.gz"
# label="gencodev45_k15"
# trait="overallBC"
# bin_index=1

####################################################################################
# begin
####################################################################################

set.seed(123)

print(paste("Running bin", bin_index))

print(label)
sumstats = fread(gwas_file)
sumstats <-  sumstats %>%
  mutate(SNP = paste(CHR, BP, A1, A2, sep = ":"))
sumstats <-  sumstats %>%
mutate(SNP2 = paste(CHR, BP, A2, A1, sep = ":"))

# Set up output_overlap_ss directory
dir.create(paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/', label, '_associations'), recursive = TRUE)
outFile_isotwas = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/', label, '_associations/', trait, '_isoTWAS.tsv')

# Identify which genes still need to be done
fff = list.files(isotwas_folder)
toDoGenes = sapply(strsplit(fff, '_isoTWAS'), function(x) x[1])

# Divide into 100 bins and select the one for this run
n_bins <- 100
bins <- split(toDoGenes, cut(seq_along(toDoGenes), n_bins, labels = FALSE))
toDoGenes <- bins[[bin_index]]

# Remove genes already completed
if (file.exists(outFile_isotwas)) {
  iii = fread(outFile_isotwas)
  toDoGenes = toDoGenes[!toDoGenes %in% iii$Gene]
}

# Sanity check
print(paste("Genes in bin", bin_index, ":", length(toDoGenes)))

doneFile = paste0('/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged/doneFile/',label,'_',trait,'_done_dev.tsv')

#file.remove(doneFile)

if (!exists('iii')){
  
  fwrite(data.frame(Trait = trait,
                    Gene = 'test',
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}

if (!file.exists(doneFile)){
  
  fwrite(data.frame(Trait = trait,
                    Gene = unique(c(iii$Gene)),
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}

runAssociation = function(gene){
  
  print(gene)
  
  if (nrow(subset(vroom::vroom(doneFile,show_col_types=F),
                  Gene == gene & Trait == trait)) == 0){
    
    ### RUN ISOTWAS
    isotwas_model = readRDS(file.path(isotwas_folder,
                                      paste0(gene,'_isoTWAS.RDS')))
    isotwas_model$R2 = unlist(isotwas_model$R2)
    isotwas_model$rsID = isotwas_model$SNP
    isotwas_model$SNP = paste0(isotwas_model$Chromosome, ":", isotwas_model$Position, ":", isotwas_model$REF, ":", isotwas_model$ALT)
    isotwas_model$SNP <- gsub("chr", "", isotwas_model$SNP)  # Remove "chr"
    

    if (any(unlist(isotwas_model$R2) > 0.01)){
      
      LD_isotwas = readRDS(file.path(ld_folder,
                                     paste0(gene,'_LDMatrix.RDS')))
      current_names <- rownames(LD_isotwas)
      
      # Create named vector for mapping
      rs_to_snp <- setNames(isotwas_model$SNP, isotwas_model$rsID)
      
      # Map current names to new SNP_IDs
      new_names <- rs_to_snp[current_names]
      
      # Check for missing values
      if (any(is.na(new_names))) {
        warning("Some rsIDs could not be mapped to SNP_IDs. Check your mapping data.frame.")
      }
      
      # Assign new names
      rownames(LD_isotwas) <- new_names
      colnames(LD_isotwas) <- new_names
      
      # Copy so we don't overwrite original
      sumstats2 <- copy(sumstats)

            # Make sure both are data.tables
      setDT(sumstats2)
      setDT(isotwas_model)


      # --- Step 1: mark matches in SNP and SNP2 ---
      setkey(sumstats2, SNP)
      setkey(isotwas_model, SNP)

      # matches in SNP
      matched_snp <- sumstats2[SNP %in% isotwas_model$SNP]

      # matches in SNP2
      matched_snp2 <- sumstats2[SNP2 %in% isotwas_model$SNP]

      # --- Step 2: flip alleles + BETA for SNP2 matches ---
      # Swap alleles back so they align with isotwas_model
      matched_snp2[, c("A1","A2") := .(A2,A1)]
      matched_snp2[, BETA := -BETA]

      # --- Step 3: combine ---
      sumstats.cur <- rbind(matched_snp, matched_snp2, fill=TRUE)

      # --- Step 4: keep only unique SNPs that are in isotwas_model ---
      sumstats.cur <- sumstats.cur[SNP %in% isotwas_model$SNP | SNP2 %in% isotwas_model$SNP]
      sumstats.cur[, Z := BETA / SE]

      tot = isotwas_model[,c('SNP','Chromosome','Position')]
      tot = tot[!duplicated(tot$SNP),]
      sumstats.cur = merge(sumstats.cur,tot,by = 'SNP')
      sumstats.cur$BETA = sumstats.cur$Z
      sumstats.cur$SE = 1
      
      isotwas_model = subset(isotwas_model,
                             SNP %in% sumstats.cur$SNP)
      isotwas_model$A1 = isotwas_model$REF
      isotwas_model$A2 = isotwas_model$ALT
      isotwas_model$Transcript = isotwas_model$Feature
      print(max(isotwas_model$R2))
      
      
      isotwas_model = subset(isotwas_model,R2 >= 0.01)
      for (tx in unique(isotwas_model$Transcript)){
        
        if (nrow(isotwas_model) > 0){
          #print(tx)
        
          tx_df_isotwas = burdenTest(mod = subset(isotwas_model,
                                                  Transcript == tx),
                                     ld = LD_isotwas,
                                     gene = gene,
                                     sumStats = sumstats.cur,
                                     chr = 'Chromosome',
                                     pos = 'Position',
                                     a1 = 'A2',
                                     a2 = 'A1',
                                     Z = 'Z',
                                     beta = 'BETA',
                                     se = 'SE',
                                     R2cutoff = 0.01,
                                     alpha = 1e-3,
                                     nperms = 1e3,
                                     usePos = F)

          if (class(tx_df_isotwas) == 'data.frame'){
            
            tx_df_isotwas$R2 = subset(isotwas_model,
                                      Transcript == tx)$R2[1]
            colnames(tx_df_isotwas) = c('Gene','Transcript',
                                        'Z','P','permute.P','topSNP',
                                        'topSNP.P','R2')
            
            fwrite(tx_df_isotwas,outFile_isotwas,
                   append = T, sep = '\t',
                   quote = F, row.names=F)
            
          }
        }
        
      }
      
    }
    
    
    
    fwrite(data.frame(Trait = trait,
                      Gene = gene,
                      Index = 1),
           doneFile,
           sep='\t',
           row.names=F,
           quote=F,append=T)
  }
}


# Run association function
require(pbmcapply)
tc_assoc = function(s){
  tryCatch(runAssociation(s),
           error = function(e) {
             print(paste0('error with ',s))
           }
  )
}

# Only apply to the selected bin
lapply(toDoGenes,
       tc_assoc)
