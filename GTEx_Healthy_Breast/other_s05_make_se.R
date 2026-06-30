library(tximeta)
library(dplyr)
library(edgeR)
library(data.table)

fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/gencodev45")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/gencodev45/",
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

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/gencodev45/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_gencodev45.RData")    
}





### veigaNorm


fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/veigaNorm")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/veigaNorm/",
    sample_list$sample,"/quant.sf")

if(sum(sample_list$sample==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","sample")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_normal_filtered/ESPRESSO_normal_filtered.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/ESPRESSO_normal_filtered/index/ESPRESSO_normal_filtered",
                    source="Filtered ESPRESSO assembly of nromal samples from Veiga et al.",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/ESPRESSO_normal_filtered/index/gentrome.fa",
                    gtf=gtf,
                    write=FALSE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T)

    se.g <- summarizeToGene(se)

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/veigaNorm/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_Veiga_ESPRESSO_normal_filtered.RData")    
}




## combined viegaNorm + gencodev45

fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/combined_veigaNorm_gencodev45")

sample_list <- data.frame(fread("/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"))

rownames(sample_list) <- sample_list$sample
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/combined_veigaNorm_gencodev45/",
    sample_list$sample,"/quant.sf")

if(sum(sample_list$sample==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","sample")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_normal_filtered.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_normal_filtered/index",
                    source="Filtered ESPRESSO assembly of normal samples from Veiga et al. combined with gencodev45",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_normal_filtered/gentrome.fa",
                    gtf=gtf,
                    write=FALSE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T)

    se.g <- summarizeToGene(se)

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/combined_veigaNorm_gencodev45/",
    sample_list$sample)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/tximeta/GTEx_breast_mammary_combined_veigaNorm_gencodev45.RData")    
}

