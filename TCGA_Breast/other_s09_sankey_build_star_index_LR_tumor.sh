#!/bin/bash
#BSUB -J "build_star_indx_LR_tumor"
#BSUB -o /rsrch5/scratch/epi/sthead/sankey/logs/build_star_index_LR_tumor_%J_%I.out
#BSUB -e /rsrch5/scratch/epi/sthead/sankey/logs/build_star_index_LR_tumor_%J_%I.out
#BSUB -q medium
#BSUB -W 4:00
#BSUB -n 2
#BSUB -M 50
#BSUB -R rusage[mem=50]

GENOME_FA=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta
#GTF_GENCODE=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.gtf
GTF_LR=/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/ESPRESSO_tumor_filtered/ESPRESSO_tumor_filtered.gtf

OUT_DIR=/rsrch5/scratch/epi/sthead/sankey/star_index_LR_tumor

mkdir -p ${OUT_DIR}

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1
module add star

# Build STAR index (GENCODE)
STAR --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir $OUT_DIR \
  --genomeFastaFiles $GENOME_FA \
  --sjdbGTFfile $GTF_LR \
  --sjdbOverhang 49
