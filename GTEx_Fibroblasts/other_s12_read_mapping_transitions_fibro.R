library(dplyr)
library(plyr)
library(ggplot2)
library(ggalluvial)
library(data.table)
library(readr)
library(purrr)
library(rtracklayer)

####################################################################################
# parse and set arguments
####################################################################################

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
sample_ind <- as.numeric(args[1])
gene <- as.character(args[2])

# sample_ind <- 3
# gene <- "SMC2"

################################################################################
# Read Mapping Transition Analysis
################################################################################

sampledat <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"))

sample_id <- sampledat$sample[sample_ind]

data_dir <- "/rsrch5/scratch/epi/sthead/sankey/tsv" # where my .tsv files were stored

# Requires .tsv files in the format: read_id transcript_id
# Paired, for each sample, for 2 annotations (e.g. gencode and assembly)
# Generate these from BAM files with:
## samtools view input.bam | awk -F'\t' '{print $1 "\t" $3}' > read_to_transcript.tsv

## Parallel processing function for read mapping transitions
process_sample <- function(sample_id) {
  file_A <- file.path(data_dir, paste0(sample_id, "_fibro_gencode.tsv.gz"))   # gencode = Assembly A
  file_B <- file.path(data_dir, paste0(sample_id, "_fibro_LR.tsv.gz"))  # assembly = Assembly B
  
  if (!file.exists(file_A) || !file.exists(file_B)) {
    message("Skipping missing files for: ", sample_id)
    return(NULL)
  }
  
  message("Processing sample: ", sample_id)
  
  message("Importing reads_A for: ", sample_id)
  reads_A <- as.data.table(read_tsv(file_A, col_names = c("read", "transcript_A"), show_col_types = FALSE))
  
  message("Importing reads_B for: ", sample_id)
  reads_B <- as.data.table(read_tsv(file_B, col_names = c("read", "transcript_B"), show_col_types = FALSE))
  
  message("Processing transitions for: ", sample_id)

  setkey(reads_A, read)
  setkey(reads_B, read)

  reads_A_collapsed <- reads_A[
    , .(transcript_A = fifelse(
          transcript_A %in% target_A, 
          transcript_A,                     # keep the transcript name
          fifelse(transcript_A == "*", "unmapped", "Other_A")  # else label
        ),
        weight = 1 / .N                      # fractional weight per read
    ), 
    by = read
  ]

  # Collapse by read and transcript label, summing the fractional weights
  reads_A_collapsed_summed <- reads_A_collapsed[
    , .(weight = sum(weight)),
    by = .(read, transcript_A)
  ]

  rm(reads_A, reads_A_collapsed)


  reads_B_collapsed <- reads_B[
    , .(transcript_B = fifelse(
          transcript_B %in% target_B, 
          transcript_B,                     # keep the transcript name
          fifelse(transcript_B == "*", "unmapped", "Other_B")  # else label
        ),
        weight = 1 / .N                      # fractional weight per read
    ), 
    by = read
  ]

  # Collapse by read and transcript label, summing the fractional weights
  reads_B_collapsed_summed <- reads_B_collapsed[
    , .(weight = sum(weight)),
    by = .(read, transcript_B)
  ]


  # ensure column names match Sankey expectations
  setnames(reads_A_collapsed_summed, "transcript_A", "cat_A")
  setnames(reads_B_collapsed_summed, "transcript_B", "cat_B")

  # join by read
  combined <- reads_A_collapsed_summed[reads_B_collapsed_summed, on = "read", nomatch = 0]

  # Sum the weights for each transition
  sankey_dt <- combined[, .(Freq = sum(weight * i.weight)), by = .(cat_A, cat_B)]

  return(sankey_dt)
}

################################################################################
# code for a specific gene
################################################################################

gtf_lr <- data.frame(import("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered_cleaned.gtf"))
gtf_lr <- gtf_lr[gtf_lr$type=="transcript",]
gtf_gencode <- data.frame(import("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"))
gtf_gencode <- gtf_gencode[gtf_gencode$type=="transcript",]

gene_abbrev <- gene # gene name
gene <- gtf_gencode$gene_id[gtf_gencode$gene_name==gene_abbrev][1] # ensembl id
#gene_abbrev <- gtf_gencode$gene_name[gtf_gencode$gene_id==gene][1]
## Define gene transcripts from GENCODE for transition analysis, these will be on the left-hand side of the plot
target_A <- gtf_gencode$transcript_id[gtf_gencode$gene_id==gene]

# Set up transcript categories for assembly
classif <- data.frame(fread("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered_classification.txt"))
transcript_categories <- classif[, c("isoform", "structural_category")]
names(transcript_categories) <- c("transcript_B", "category_B")

# Select assembly transcripts for right-hand side of plot
target_B <- classif[classif$associated_gene == gene, "isoform"]

transcript_categories$category_B[!transcript_categories$transcript_B %in% target_B] <- paste0("non-",gene_abbrev)
transcript_categories$category_B[transcript_categories$category_B=="full-splice_match"] <- transcript_categories$transcript_B[transcript_categories$category_B=="full-splice_match"]
#transcript_categories[transcript_categories$transcript_B == "ENST00000304858.7", "category_B"] <- "ENST00000304858.7"
dt_T <- as.data.table(transcript_categories)
setkey(dt_T, transcript_B)

## Process sample
flow_counts_all <- process_sample(sample_id)

write.table(flow_counts_all,file=paste0("/rsrch5/scratch/epi/sthead/sankey/flow_counts/gtex_fibro_",gene_abbrev,"_",sample_id,".txt"),
  sep="\t", quote=F, row.names=F,col.names=T)

# flow_counts <- flow_counts_all[-which(flow_counts_all$cat_A=="Other_A" & flow_counts_all$cat_B=="Other_B"),]
# flow_counts <- flow_counts[-which(flow_counts$cat_A=="unmapped" & flow_counts$cat_B=="unmapped"),]
# flow_counts <- flow_counts[-which(flow_counts$cat_A=="unmapped" & flow_counts$cat_B=="Other_B"),]
# flow_counts <- flow_counts[-which(flow_counts$cat_A=="Other_A" & flow_counts$cat_B=="unmapped"),]
# flow_counts <- merge(flow_counts,transcript_categories,by.x="cat_B",by.y="transcript_B",all.x=T,all.y=F)

# ### Aggregate transitions across samples
# flow_counts_combined.HSPA4 <- flow_counts_all %>%
#   group_by(bucket_A, category_B) %>%
#   summarise(Freq = sum(Freq), .groups = "drop")
# #flow_counts_combined.HSPA4 <- flow_counts_all

# ### Revalue factor levels for visualization
# flow_counts_combined.HSPA4$bucket_A <- revalue(flow_counts_combined.HSPA4$bucket_A, 
#                                               c("Other_A" = paste0("non-HSPA4\n(N = ",
#                                                                    252930 - length(target_A), ")")))
# # need to check the value used in N=.. in above line 252930 is the number f transcripts in gencodev45

# flow_counts_combined.HSPA4$category_B <- revalue(flow_counts_combined.HSPA4$category_B, 
#                                                 c("non-HSPA4" = paste0("non-HSPA4\n(N = ",
#                                                                       length(transcript_categories[transcript_categories$category_B == "non-HSPA4", "transcript_B"]), ")"),
#                                                   "full-splice_match" = paste0("full-splice_match\n(N = ",
#                                                                  length(transcript_categories[transcript_categories$category_B == "full-splice_match", "transcript_B"]), ")"),
#                                                   "incomplete-splice_match" = paste0("incomplete-splice_match\n(N = ",
#                                                                  length(transcript_categories[transcript_categories$category_B == "incomplete-splice_match", "transcript_B"]), ")"),
#                                                   "novel_in_catalog" = paste0("novel_in_catalog\n(N = ",
#                                                                  length(transcript_categories[transcript_categories$category_B == "novel_in_catalog", "transcript_B"]), ")"),
#                                                   "novel_not_in_catalog" = paste0("novel_not_in_catalog\n(N = ",
#                                                                  length(transcript_categories[transcript_categories$category_B == "novel_not_in_catalog", "transcript_B"]), ")")))

# flow_counts_combined.HSPA4 <- data.frame(flow_counts_combined.HSPA4)

# # Prepare data for visualization (remove unmapped and low-frequency transitions)
# flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4
# rm <- intersect(grep("non-",flow_counts_combined.HSPA4.plot$bucket_A),grep("non-",flow_counts_combined.HSPA4.plot$category_B))
# rm <- c(rm, intersect(grep("unmapped",flow_counts_combined.HSPA4.plot$bucket_A),grep("non-",flow_counts_combined.HSPA4.plot$category_B)))
# flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4.plot[-rm, ]
# flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4.plot[-grep("unmapped", flow_counts_combined.HSPA4.plot$category_B), ]

# ## Figure 3D: HSPA4 read mapping transitions
# flow_counts <- flow_counts[flow_counts$Freq < 1,]
# fig <- ggplot(flow_counts, 
#                 aes(axis1 = cat_A, axis2 = cat_B, y = sqrt(Freq))) +
#   geom_alluvium(aes(fill = cat_B), width = 1/4) +
#   #scale_fill_manual(values = c("red", myPalette[c(1, 3, 4)], "grey30")) +
#   geom_stratum(width = 1/4, fill = "gray95", color = "black") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
#             size = 3, lineheight = 0.7) +
#   scale_x_discrete(limits = c("gencode", "lr-assembly"), expand = c(.1, .1)) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(title = paste0(gene_abbrev, " read mapping transitions"), 
#        y = expression(sqrt("Weighted read sum")), x = NULL)

# ggsave(paste0(gene_abbrev,"_sankey.png"), fig, width = 10, height = 10, units = "in", dpi = 300, bg = 'white')
