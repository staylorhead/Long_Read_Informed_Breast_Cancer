library(dplyr)
library(plyr)
library(ggplot2)
library(ggalluvial)
library(data.table)
library(readr)
library(purrr)

################################################################################
# Read Mapping Transition Analysis
################################################################################

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
  
  combined <- data.table(transcript_A = reads_A$transcript_A,
                         transcript_B = reads_B$transcript_B)
  
  rm(reads_A, reads_B)
  
  combined[, bucket_A := fifelse(
    transcript_A == "*", "unmapped",
    fifelse(transcript_A %in% target_A, transcript_A, "Other_A")
  )]
  
  setkey(combined, transcript_B)
  final <- dt_T[combined]
  rm(combined)
  
  final[is.na(category_B), category_B := "unmapped"]
  
  reads_summary <- data.frame(table(as.data.frame(final[, .(bucket_A, category_B)])))
  
  return(reads_summary)
}

# Here I selected a random set of 5 samples 

sampledat <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"))

sample_ids <- sampledat$sample[1:5]

data_dir <- "/rsrch5/scratch/epi/sthead/sankey/tsv" # where my .tsv files were stored

################################################################################
# Figure 3D: HSPA4 Read Mapping Transitions
################################################################################

## Define HSPA4 transcripts from GENCODE for transition analysis, these will be on the left-hand side of the plot
target_A <- c('ENST00000304858.7','ENST00000504328.1','ENST00000514825.1')

# Set up transcript categories for assembly
classif <- data.frame(fread("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered_classification.txt"))



transcript_categories <- classif[, c("isoform", "structural_category")]


transcript_categories$structural_category <- as.character(transcript_categories$structural_category)
names(transcript_categories) <- c("transcript_B", "category_B")

# Select assembly transcripts for right-hand side of plot
target_B <- classif[classif$associated_gene == "ENSG00000170606.16", "isoform"]
# Here I separated the transcripts by structural category but highlighted one isoform of interest
# Also noted non-HSPA4 isoforms
transcript_categories[!transcript_categories$transcript_B %in% target_B, "category_B"] <- "non-HSPA4"
transcript_categories[transcript_categories$transcript_B == "ENST00000304858.7", "category_B"] <- "ENST00000304858.7"
dt_T <- as.data.table(transcript_categories)
setkey(dt_T, transcript_B)

## Process GDM = 0 samples
flow_counts_all <- map_dfr(sample_ids, process_sample)

### Aggregate transitions across samples
flow_counts_combined.HSPA4 <- flow_counts_all %>%
  group_by(bucket_A, category_B) %>%
  summarise(Freq = sum(Freq), .groups = "drop")
#flow_counts_combined.HSPA4 <- flow_counts_all

### Revalue factor levels for visualization
flow_counts_combined.HSPA4$bucket_A <- revalue(flow_counts_combined.HSPA4$bucket_A, 
                                              c("Other_A" = paste0("non-HSPA4\n(N = ",
                                                                   252930 - length(target_A), ")")))
# need to check the value used in N=.. in above line 252930 is the number f transcripts in gencodev45

flow_counts_combined.HSPA4$category_B <- revalue(flow_counts_combined.HSPA4$category_B, 
                                                c("non-HSPA4" = paste0("non-HSPA4\n(N = ",
                                                                      length(transcript_categories[transcript_categories$category_B == "non-HSPA4", "transcript_B"]), ")"),
                                                  "full-splice_match" = paste0("full-splice_match\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "full-splice_match", "transcript_B"]), ")"),
                                                  "incomplete-splice_match" = paste0("incomplete-splice_match\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "incomplete-splice_match", "transcript_B"]), ")"),
                                                  "novel_in_catalog" = paste0("novel_in_catalog\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "novel_in_catalog", "transcript_B"]), ")"),
                                                  "novel_not_in_catalog" = paste0("novel_not_in_catalog\n(N = ",
                                                                 length(transcript_categories[transcript_categories$category_B == "novel_not_in_catalog", "transcript_B"]), ")")))

flow_counts_combined.HSPA4 <- data.frame(flow_counts_combined.HSPA4)

# Prepare data for visualization (remove unmapped and low-frequency transitions)
flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4
rm <- intersect(grep("non-",flow_counts_combined.HSPA4.plot$bucket_A),grep("non-",flow_counts_combined.HSPA4.plot$category_B))
rm <- c(rm, intersect(grep("unmapped",flow_counts_combined.HSPA4.plot$bucket_A),grep("non-",flow_counts_combined.HSPA4.plot$category_B)))
flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4.plot[-rm, ]
flow_counts_combined.HSPA4.plot <- flow_counts_combined.HSPA4.plot[-grep("unmapped", flow_counts_combined.HSPA4.plot$category_B), ]

## Figure 3D: HSPA4 read mapping transitions
fig <- ggplot(flow_counts_combined.HSPA4.plot, 
                aes(axis1 = bucket_A, axis2 = category_B, y = sqrt(Freq))) +
  geom_alluvium(aes(fill = category_B), width = 1/4) +
  #scale_fill_manual(values = c("red", myPalette[c(1, 3, 4)], "grey30")) +
  geom_stratum(width = 1/4, fill = "gray95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, lineheight = 0.7) +
  scale_x_discrete(limits = c("gencode", "lr-assembly"), expand = c(.1, .1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "HSPA4 Read Mapping Transitions", 
       y = expression(sqrt("Read Count")), x = NULL)

ggsave("HSPA4_sankey.png", fig, width = 7.74, height = 6.17, units = "in", dpi = 300, bg = 'white')