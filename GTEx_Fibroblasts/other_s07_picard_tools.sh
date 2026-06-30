#!/bin/sh
#BSUB -J picard[1-504]
#BSUB -o /rsrch5/home/epi/sthead/isoqtl_lr_breast/logs/other/picard_fibro_%J_%I.out
#BSUB -e /rsrch5/home/epi/sthead/isoqtl_lr_breast/logs/other/picard_fibro_%J_%I.err
#BSUB -W 3:00
#BSUB -n 5
#BSUB -M 20
#BSUB -R rusage[mem=20]
#BSUB -q short

DIR_BAM=/rsrch5/scratch/epi/sthead/GTEx_v8/bam
SAMPLES_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_fibroblast_bam_files_with_sample_attrib.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)

OUT_DIR="/rsrch5/scratch/epi/sthead/GTEx_v8/picard/${SAMPLE}"
FA_FILE="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#TMP_DIR="/rsrch5/scratch/epi/sthead/GTEx_v8/tmp_bam/${SAMPLE}"
REFFLAT='/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.nochr.refflat'

mkdir -p ${OUT_DIR}
#mkdir -p ${TMP_DIR}

module load picard

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate picard-2.27.4

picard CollectGcBiasMetrics \
    I=${DIR_BAM}/${SAMPLE}.sortedbyCoord_Aligned.sortedByCoord.out.bam \
    O=${OUT_DIR}/gcBiasMetrics.txt \
    S=${OUT_DIR}/gcSummaryMetrics.txt \
    CHART=${OUT_DIR}/gcBiasMetrics.pdf \
    R=${FA_FILE}

picard CollectRnaSeqMetrics \
    INPUT=${DIR_BAM}/${SAMPLE}.sortedbyCoord_Aligned.sortedByCoord.out.bam \
    OUTPUT=${OUT_DIR}/txomeMetrics.txt \
    REF_FLAT=${REFFLAT} \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND    

picard CollectAlignmentSummaryMetrics \
    INPUT=${DIR_BAM}/${SAMPLE}.sortedbyCoord_Aligned.sortedByCoord.out.bam \
    OUTPUT=${OUT_DIR}/genomeMetrics.txt \
    REFERENCE_SEQUENCE=${FA_FILE}

