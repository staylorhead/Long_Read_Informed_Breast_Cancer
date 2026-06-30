library(data.table)
library(stringr)
bams <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_bams.txt",header=F)
bams_full <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_bams_full_filename.txt",header=F)

colnames(bams_full) <- "filename"
# Extract just the sample name (everything before the first extension)
bams_full$sample <- sub("\\..*$", "", bams_full$filename)

# Extract the extension (everything after the first dot)
bams_full$extension <- sub("^.*?\\.", "", bams_full$filename)
sum(bams_full$sample==bams$V1) # good

sample_attribs <- fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                        header=T)
sum(bams_full$sample %in% sample_attribs$SAMPID)
length(unique(bams_full$sample))

comb <- merge(x=bams_full,y=sample_attribs,by.x="sample",by.y="SAMPID",all.x=T,all.y=F)
table(comb$SMTSD) # 525 cultured fibroblasts
fibro <- data.frame(comb[comb$SMTSD=="Cells - Cultured fibroblasts",])
fibro <- fibro[,c(2,1,3:ncol(fibro))]
fibro <- fibro[fibro$SMAFRZE=="RNASEQ",]
fibro$SUBJID <- sapply(strsplit(fibro$sample, "-"), function(x) paste(x[1:2], collapse = "-"))
length(unique(fibro$SUBJID)) # 504

write.table(fibro,file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt",
            sep="\t",row.names = F,col.names = T,quote=F)
