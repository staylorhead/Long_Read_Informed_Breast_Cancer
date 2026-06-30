# Packages
.libPaths(c( "/home/stbresnahan/R/ubuntu/4.3.1" , .libPaths()))
library(plyr)
library(dplyr)
library(seqinr)
library(rtracklayer)
library(tidyr)
library(doParallel)
library(foreach)


# Read in tumor assembly annotation as data.frame
gtf_tumor <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_tumor/ESPRESSO_tumor_corrected.gtf"
gtf_tumor <- data.frame(readGFF(gtf_tumor))

# Read in tumor assembly isoform classification table and format
classif <- read.table("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_tumor/ESPRESSO_tumor_classification.txt",header=T)
names(classif)[1] <- "transcript_id"
classif$structural_category <- factor(as.character(classif$structural_category),
                                      levels=c("full-splice_match",
                                               "incomplete-splice_match",
                                               "novel_in_catalog",
                                               "novel_not_in_catalog",
                                               "genic",
                                               "antisense",
                                               "fusion",
                                               "intergenic",
                                               "genic_intron"))
classif$structural_category <- revalue(classif$structural_category, 
                                       c("full-splice_match"="FSM",
                                         "incomplete-splice_match"="ISM",
                                         "novel_in_catalog"="NIC",
                                         "novel_not_in_catalog"="NNC",
                                         "genic"="Genic Genomic",
                                         "antisense"="Antisense",
                                         "fusion"="Fusion",
                                         "intergenic"="Intergenic",
                                         "genic_intron"="Genic Intron"))

# Read in filtered tumor assembly classification table and format
classif_filtered <- read.table("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_tumor_filtered/ESPRESSO_tumor_filtered_classification.txt",header=T)
names(classif_filtered)[1] <- "transcript_id"

# Set up tx2g, setting gene_id to classif$associated_gene if missing
annot.tumor <- gtf_tumor
annot.tumor$score <- NULL # remove
annot.tumor$phase <- NULL # remove
novelTx.tumor <- unique(annot.tumor[annot.tumor$gene_id=="NA","transcript_id"])
novelTx.ascg <- classif[classif$transcript_id%in%novelTx.tumor,c("transcript_id","associated_gene")]
novelTx.ascg <- novelTx.ascg[grep("ENSG",novelTx.ascg$associated_gene),]
annot.tumor <- left_join(annot.tumor,novelTx.ascg,by="transcript_id")
annot.tumor[is.na(annot.tumor$associated_gene),"associated_gene"] <- "NA"
annot.tumor[annot.tumor$gene_id=="NA","gene_id"] <- annot.tumor[annot.tumor$gene_id=="NA","associated_gene"]
annot.tumor$associated_gene <- NULL

# Create gene_ids for novel genes
novelTx.tumor.unannot <- setdiff(novelTx.tumor,novelTx.ascg$transcript_id)
annot.tumor[annot.tumor$transcript_id%in%novelTx.tumor.unannot,"gene_id"] <- sapply(annot.tumor[annot.tumor$transcript_id%in%novelTx.tumor.unannot,"transcript_id"],
                                                                                    function(x) paste(strsplit(x,":")[[1]][1:3],collapse=":"))
tx2g.tumor <- annot.tumor[,c("transcript_id","gene_id")]
tx2g.tumor <- tx2g.tumor[!duplicated(tx2g.tumor),]

# left join gtf_tumor with tx2g.tumor so that gene_ids match
gtf_tumor$score <- "."
gtf_tumor$phase <- "."
gtf_tumor$gene_id <- NULL
gtf_tumor <- left_join(gtf_tumor,tx2g.tumor,by="transcript_id")
gtf_tumor <- gtf_tumor[,c(1:9,11,10)]
gtf_tumor$info <- NA

# Subset to the high-confidence isoforms
gtf_tumor <- gtf_tumor[gtf_tumor$transcript_id%in%classif_filtered$transcript_id,]

# Format back to gtf
gtf_tumor[gtf_tumor$type=="transcript","info"] <- paste0("transcript_id \"",
                                                         gtf_tumor[gtf_tumor$type=="transcript","transcript_id"],
                                                         "\"; gene_id \"",
                                                         gtf_tumor[gtf_tumor$type=="transcript","gene_id"],"\"")
gtf_tumor[gtf_tumor$type=="exon","info"] <- paste0("transcript_id \"",
                                                   gtf_tumor[gtf_tumor$type=="exon","transcript_id"],
                                                   "\"; gene_id \"",
                                                   gtf_tumor[gtf_tumor$type=="exon","gene_id"],"\";")
gtf_tumor$exon_number <- NULL
gtf_tumor <- gtf_tumor[,c(10,9,1:8,11)]

# Read in GENCODE v45 annotation as data.frame and wrangle to same format as gtf_tumor
## Note: tthis file is the same as /rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf
### subset to chr 1-22 X, Y, M
gtf_gencode <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/gencode_v45_dedup.gtf"
gtf_gencode <- data.frame(readGFF(gtf_gencode))
gtf_gencode <- gtf_gencode[!is.na(gtf_gencode$transcript_id),] # Remove lines with missing transcript ids
gtf_gencode$type <- as.character(gtf_gencode$type)
gtf_gencode <- gtf_gencode[gtf_gencode$type%in%c("transcript","exon"),c("gene_id","transcript_id","seqid","source","type",
                                                                               "start","end","score","strand","phase")]
gtf_gencode$score <- "."
gtf_gencode$phase <- "."

gtf_gencode$info <- NA
gtf_gencode[gtf_gencode$type=="transcript","info"] <- paste0("transcript_id \"",
                                                             gtf_gencode[gtf_gencode$type=="transcript","transcript_id"],
                                                             "\"; gene_id \"",
                                                             gtf_gencode[gtf_gencode$type=="transcript","gene_id"],"\"")
gtf_gencode[gtf_gencode$type=="exon","info"] <- paste0("transcript_id \"",
                                                       gtf_gencode[gtf_gencode$type=="exon","transcript_id"],
                                                       "\"; gene_id \"",
                                                       gtf_gencode[gtf_gencode$type=="exon","gene_id"],"\";")
gtf_gencode <- gtf_gencode[,names(gtf_tumor)]

# Get genes and transcripts in gtf_tumor that are not in gtf_gencode
novel_genes <- setdiff(gtf_tumor$gene_id,gtf_gencode$gene_id)
novel_transcripts <- setdiff(gtf_tumor$transcript_id,gtf_gencode$transcript_id)
novel_transcripts <- novel_transcripts[!is.na(novel_transcripts)]

# Combine gtf_gencode with novel transcripts and genes in gtf_tumor
gtf_combined <- rbind(gtf_gencode,gtf_tumor[gtf_tumor$gene_id%in%novel_genes,])
gtf_combined <- rbind(gtf_combined,gtf_tumor[gtf_tumor$transcript_id%in%novel_transcripts,])
gtf_combined <- gtf_combined[!duplicated(gtf_combined),]
gtf_combined$source <- as.character(gtf_combined$source)

# Set source value
gtf_combined[gtf_combined$source%in%c("HAVANA","ENSEMBL"),"source"] <- "GENCODEv45"
gtf_combined[gtf_combined$source%in%c("novel","novel_isoform"),"source"] <- "tumor_assembly"

# Sort so that transcripts appear after genes and exons are in order
gtf_sorted <- gtf_combined
gtf_sorted$original_index <- 1:nrow(gtf_sorted)
gtf_sorted$group_key <- gtf_sorted$transcript_id
gtf_sorted$type_order <- ifelse(gtf_sorted$type == "transcript", 1, 2)
gtf_sorted <- gtf_sorted[order(gtf_sorted$seqid,gtf_sorted$start, gtf_sorted$group_key, gtf_sorted$type_order), ]
gtf_sorted$original_index <- NULL
gtf_sorted$group_key <- NULL
gtf_sorted$type_order <- NULL
gtf_sorted$gene_id <- NULL
gtf_sorted$transcript_id <- NULL
gtf_combined <- gtf_sorted

# Save combined gtf
file.combined <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_tumor_filtered.gtf"
write.table(gtf_combined_out,file.combined,
            row.names=F,col.names=F,quote=F,sep="\t")



# Add top level (gene) entries if needed
gtf_combined_in <- data.frame(import(file.combined))
gtf_combined_in$score <- "."
gtf_combined_in$phase <- "."

# Function to add gene entries (coordinates == longest transcript)
add_top_entries_parallel <- function(df, num_cores = detectCores() - 1) {
  # Set up parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Get unique gene IDs
  unique_genes <- unique(df$gene_id)
  
  # Process genes in parallel
  results <- foreach(gene = unique_genes, .combine = rbind, .packages = "base") %dopar% {
    # Get all rows for this gene
    gene_rows <- df[df$gene_id == gene, ]
    
    # Find the FIRST entry (first occurrence of this gene_id)
    first_entry <- gene_rows[1, ]
    
    # Calculate start and end for TOP entry
    min_start <- min(gene_rows$start)
    max_end <- max(gene_rows$end)
    
    # Create TOP entry
    top_entry <- first_entry
    if(any(gene_rows$source%in%c("GENCODEv45"))){
      top_entry$source <- "GENCODEv45"
      }else{
      top_entry$source <- "tumor_assembly"
    }
    top_entry$type <- "gene"
    top_entry$start <- min_start
    top_entry$end <- max_end
    top_entry$width <- abs(max_end-min_start)+1
    top_entry$strand <- first_entry$strand
    top_entry$transcript_id <- NA
    
    # Return TOP entry followed by all gene rows
    rbind(top_entry, gene_rows)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Reset row names
  rownames(results) <- NULL
  
  return(results)
}

# Add top level entires
gtf_combined_out <- add_top_entries_parallel(gtf_combined_in)

# Convert back to gtf
gtf_combined_out$info <- NA
gtf_combined_out[gtf_combined_out$type=="transcript","info"] <- paste0("transcript_id \"",
                                                                       gtf_combined_out[gtf_combined_out$type=="transcript","transcript_id"],
                                                                       "\"; gene_id \"",
                                                                       gtf_combined_out[gtf_combined_out$type=="transcript","gene_id"],"\"")
gtf_combined_out[gtf_combined_out$type=="exon","info"] <- paste0("transcript_id \"",
                                                                 gtf_combined_out[gtf_combined_out$type=="exon","transcript_id"],
                                                                 "\"; gene_id \"",
                                                                 gtf_combined_out[gtf_combined_out$type=="exon","gene_id"],"\";")

# Save combined gtf
write.table(gtf_combined_out,"/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_tumor_filtered_WithGeneIDs.gtf",
            row.names=F,col.names=F,quote=F,sep="\t")

