library(tximeta)
library(dplyr)
library(edgeR)
library(data.table)

# gencode

fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_gencodev45")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_gencodev45/",
    sample_list$sample,"/quant.sf")

if(sum(sample_list$sample==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","sample")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45",
                    source="GENCODEv45",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
                    gtf=gtf,
                    write=FALSE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T)

    se.g <- summarizeToGene(se)

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_gencodev45/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_fibroblast_gencodev45.RData")    
}





### fibroblast LR


fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_ESPRESSO")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_ESPRESSO/",
    sample_list$sample,"/quant.sf")

if(sum(sample_list$sample==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","sample")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered_cleaned.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/orthogonal/Expression/short-read/salmon/GTEx_fibroblast/ESPRESSO_fibroblast_filtered/index",
                    source="Filtered ESPRESSO assembly of GTEx V9 fibroblast",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/orthogonal/Expression/short-read/salmon/GTEx_fibroblast/ESPRESSO_fibroblast_filtered/gentrome.fa",
                    gtf=gtf,
                    write=FALSE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T)

    se.g <- summarizeToGene(se)

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_ESPRESSO/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_fibroblast_ESPRESSO_filtered.RData")    
}




## combined fibroblast LR + gencodev45

fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_combined_gencodev45_and_ESPRESSO")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_combined_gencodev45_and_ESPRESSO/",
    sample_list$sample,"/quant.sf")

if(sum(sample_list$sample==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","sample")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/combined_gencodev45_and_ESPRESSO_fibroblast_filtered.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/orthogonal/Expression/short-read/salmon/GTEx_fibroblast/combined_gencodev45_and_ESPRESSO_fibroblast_filtered/index",
                    source="Filtered ESPRESSO assembly of fibroblast samples from GTEx V9 et al. combined with gencodev45",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/orthogonal/Expression/short-read/salmon/GTEx_fibroblast/combined_gencodev45_and_ESPRESSO_fibroblast_filtered/gentrome.fa",
                    gtf=gtf,
                    write=TRUE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T,skipMeta=T)

    #se.g <- summarizeToGene(se) # not working for combined annotation

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/fibroblast_combined_gencodev45_and_ESPRESSO/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    #save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_fibroblast_combined_ESPRESSO_gencodev45.RData")    
    save(se,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_fibroblast_combined_ESPRESSO_gencodev45.RData")    

}

