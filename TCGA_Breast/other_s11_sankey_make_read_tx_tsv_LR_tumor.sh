#!/bin/bash
#BSUB -J "make_tsv[1-5]"
#BSUB -o /rsrch5/scratch/epi/sthead/sankey/logs/make_tsv_LR_tumor_%J_%I.out
#BSUB -e /rsrch5/scratch/epi/sthead/sankey/logs/make_tsv_LR_tumor_%J_%I.out
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 20
#BSUB -R rusage[mem=20]

# LR TUMOR ALIGNED RESULTS

INPUT_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/sample_lists/sample_list_postQC_final_for_analysis_TCGA_present_all_annot_with_fastq_hash.txt"
ROWID=$LSB_JOBINDEX
HASH=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f2)
SAMPLE=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f3)

module load samtools

DIR_ALIGN="/rsrch5/scratch/epi/sthead/sankey/star_alignments_LR_tumor"
DIR_OUT="/rsrch5/scratch/epi/sthead/sankey/tsv"

mkdir -p ${DIR_OUT}

samtools view ${DIR_ALIGN}/${SAMPLE}_Aligned.toTranscriptome.out.bam \
  | awk -F'\t' '{print $1 "\t" $3}' \
  | gzip > ${DIR_OUT}/${SAMPLE}_tumor_LR.tsv.gz