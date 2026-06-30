#!/bin/sh
#BSUB -J bam_to_fastq[1-440]%25
#BSUB -o /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/bam_to_fastq_%J_%I.out
#BSUB -e /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/bam_to_fastq_%J_%I.err
#BSUB -W 3:00
#BSUB -n 5
#BSUB -M 50
#BSUB -R rusage[mem=50]
#BSUB -q short

DIR_BAM=/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam
SAMPLES_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)
THREADS=5
OUT_DIR_SORTED="/rsrch5/scratch/epi/sthead/GTEx_v8/bams_read_sorted"
OUT_DIR_FASTQ="/rsrch5/scratch/epi/sthead/GTEx_v8/fastq" # I also copied these over to /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/fastq

mkdir -p ${OUT_DIR_SORTED}
mkdir -p ${OUT_DIR_FASTQ}

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1

samtools sort -n -@ ${THREADS} -o ${OUT_DIR_SORTED}/${LIBID}_sorted.bam ${DIR_BAM}/${LIBID}

samtools fastq -@ ${THREADS} -1 ${OUT_DIR_FASTQ}/${LIBID}_R1.fastq -2 ${OUT_DIR_FASTQ}/${LIBID}_R2.fastq -0 /dev/null -s /dev/null -n ${OUT_DIR_SORTED}/${LIBID}_sorted.bam

rm ${OUT_DIR_SORTED}/${LIBID}_sorted.bam

conda deactivate

