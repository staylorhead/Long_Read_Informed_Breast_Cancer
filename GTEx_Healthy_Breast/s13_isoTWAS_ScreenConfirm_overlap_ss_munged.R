require(data.table)
require(isotwas)
library(rtracklayer)
library(dplyr)

trait <- c("HER2enriched","TN","lumA","lumBHER2neg","overallBC","lumB","PMD","DA","NDA","BCsurvival")

dirs <- list.dirs("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/GTEx/output_overlap_ss_munged")
dirs <- dirs[grep("associations",dirs)]

# read in transcript features variance
sd_dat <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/gene_and_tx_variance_in_bed.txt")
sd_dat <- sd_dat[sd_dat$var==0,]
sd_dat <- sd_dat[sd_dat$dataset=="GTEx",]

for (quan in c("gencodev45","veigaNorm","comb")) {
  
  if(quan=="gencodev45"){
          gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"
        }else if(quan=="comb"){
          gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_normal_filtered.gtf"
          }else if(quan=="veigaNorm"){
            gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_normal_filtered/ESPRESSO_normal_filtered.gtf"
          }

  gtf <- import(gtf_file)
  gtf <- data.frame(gtf)
  gtf <- gtf[gtf$type=="transcript",]
  gtf <- gtf[,c("transcript_id","gene_id","seqnames","start","end","width","strand","source")]

  sd_dat_quant <- sd_dat[sd_dat$annot==quan,]

  for (cancer in trait){
    
    print(quan)
    cor_dir <- dirs[grep(quan,dirs)]
    setwd(cor_dir)
    print(cancer)
    
    isotwas_res = vroom::vroom(list.files(pattern = paste0(cancer,"_"))[1], show_col_types = FALSE)
    isotwas_res = isotwas_res[complete.cases(isotwas_res),]
    isotwas_res <- isotwas_res[!isotwas_res$Transcript %in% sd_dat_quant$pid,]
    isotwas_res = isotwas_res[order(abs(isotwas_res$Z),decreasing = T),]
    isotwas_res$Gene <- sub("\\..*", "", isotwas_res$Gene)
    isotwas_res = isotwas_res[!duplicated(isotwas_res$Transcript) & 
                                abs(isotwas_res$Z) < Inf,]

    isotwas_res = merge(gtf,isotwas_res,by.x='transcript_id',by.y="Transcript",all.x=F,all.y=T)
    
    www = which(isotwas_res$seqnames == "chr6" &
                  isotwas_res$start < 35e6 &
                  isotwas_res$end > 27e6)
    isotwas_res = isotwas_res[-www,]

    # stopped here
    
    require(tidyverse)
    gene <- isotwas_res %>%
      group_by(Gene) %>%
      summarise(
        Chromosome = first(seqnames),   # assumes all same, else check consistency
        Start      = min(start, na.rm = TRUE),
        End        = max(end, na.rm = TRUE),
        Screen.P   = isotwas::p_screen(P),
      )
    
    alpha1=0.05
    G = nrow(gene)
    gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'fdr')
    R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
    alpha2 = (R*alpha1)/G
    isoform_new = as.data.frame(matrix(nrow = 0,
                                       ncol = ncol(isotwas_res)+2))
    colnames(isoform_new) = c(colnames(isotwas_res),'Screen.P','Confirmation.P')
    gene = gene[order(gene$Screen.P),]
    ttt = merge(isotwas_res,
                gene[,c('Gene','Screen.P',
                        'Screen.P.Adjusted')])
    isoform_new = ttt %>%
    group_by(Gene) %>%
    mutate(Confirmation.P = isotwas::p_confirm(P, alpha = alpha2)) %>%
    ungroup()

    isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                        isoform_new$Confirmation.P,
                                        1)
    isoform_new = isoform_new[!duplicated(isoform_new),]
    isoform_new$Indication = cancer
    isoform_sig = subset(isoform_new, 
                         Screen.P < alpha1 &
                           Confirmation.P < alpha2 &
                           permute.P < 0.05)
    
    fwrite(isoform_new,
           'ScreenConfirm_isoTWAS.tsv',
           append=T,sep="\t",
           row.names=F,
           quote=F)
    print(nrow(isoform_sig))
    
    if (nrow(isoform_sig) > 1){
      
      isoform_sig = isoform_sig[order(isoform_sig$seqnames,
                                      isoform_sig$start),]            
      fwrite(isoform_sig,
             'SignificantAssociations_isoTWAS_noFineMap.tsv',
             append=T,sep="\t",
             row.names=F,
             quote=F)
    }
    
  }
}  
