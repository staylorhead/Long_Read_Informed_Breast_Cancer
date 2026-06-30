library(tximeta)
library(dplyr)
library(edgeR)
library(data.table)

fnames <- list.files("/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered/run/TCGA")

sample_list <- data.frame(fread("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/fastq_hash_codes_to_TCGA_IDs.txt",header=F))
rownames(sample_list) <- sample_list$V1
sample_list <- sample_list[fnames,]

sample_list$quant_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered/run/TCGA/",
    sample_list$V1,"/quant.sf")

if(sum(sample_list$V1==fnames)==length(fnames)){
    metadata <- sample_list[,c("quant_path","V1")]
    names(metadata) <- c("files","names")

    setTximetaBFC("/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/tximeta")

    gtf="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/combined_gencodev45_and_ESPRESSO_tumor_filtered.gtf"

    makeLinkedTxome(indexDir="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered/index",
                    source="Filtered ESPRESSO assembly of tumor samples from Veiga et al. merged with GENCODEv45",
                    organism="Homo sapiens",
                    genome="GRCh38",
                    release="p14",
                    fasta="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered/index/gentrome.fa",
                    gtf=gtf,
                    write=FALSE)

    se <- tximeta(metadata,type="salmon",useHub=FALSE,skipSeqinfo=T,txOut=T,dropInfReps=T)

    se.g <- summarizeToGene(se)

    metadata$directories <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered/run/TCGA/",
    sample_list$V1)

    cs <- catchSalmon(metadata$directories)

    save(se,se.g,cs,file="/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/tximeta/combined_gencodev45_and_ESPRESSO_tumor_filtered.RData")    
}

