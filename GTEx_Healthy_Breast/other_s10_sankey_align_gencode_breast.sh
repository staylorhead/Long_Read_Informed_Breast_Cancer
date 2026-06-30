#!/bin/bash
#BSUB -J "align_gencode[1-5]"
#BSUB -o /rsrch5/scratch/epi/sthead/sankey/logs/align_star_gencode_breast_%J_%I.out
#BSUB -e /rsrch5/scratch/epi/sthead/sankey/logs/align_star_gencode_breast_%J_%I.out
#BSUB -q medium
#BSUB -W 5:00
#BSUB -n 8
#BSUB -M 80
#BSUB -R rusage[mem=80]

INPUT_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f2)

DIR_FASTQ="/rsrch5/scratch/epi/sthead/GTEx_v8/fastq"
DIR_ALIGN="/rsrch5/scratch/epi/sthead/sankey/star_alignments_gencode"
TXINDEX_STAR="/rsrch5/scratch/epi/sthead/sankey/star_index_gencode"
DIR_TEMP="/rsrch5/scratch/epi/sthead/sankey/temp"

# Create output directories
mkdir -p ${DIR_ALIGN}
mkdir -p ${DIR_TEMP}

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1
module add star

STAR --genomeDir ${TXINDEX_STAR} \
    --readFilesIn ${DIR_FASTQ}/${LIBID}_R1.fastq ${DIR_FASTQ}/${LIBID}_R2.fastq \
    --runThreadN 8 \
    --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate \
    --sjdbScore 1 --outTmpDir ${DIR_TEMP}/${SAMPLE} \
    --outFileNamePrefix ${DIR_ALIGN}/${SAMPLE}_ \
    --outBAMsortingBinsN 200 \
    --limitBAMsortRAM 80000000000 \
    --quantMode TranscriptomeSAM

# Clean up temporary and intermediate files
rm -rf ${DIR_TEMP}/${SAMPLE}

