#!/bin/sh
#BSUB -J gtex_quant[1-440]%100
#BSUB -o /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/quant_from_fastq_gencodev45_%J_%I.out
#BSUB -e /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/quant_from_fastq_gencodev45_%J_%I.err
#BSUB -W 3:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q short

DIR_FASTQ=/rsrch5/scratch/epi/sthead/GTEx_v8/fastq
DIR_OUT=/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/salmon_quantifications/gencodev45
TXINDEX_salmon=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v45.salmon_index/gencode_v45
THREADS=10

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

conda activate salmon-1.10.2 

INPUT_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$INPUT_FILE" | cut -f2)

# extension=Aligned.sortedByCoord.out.patched.md.bam 
# or Aligned.sortedByCoord.out.patched.md.bam,
# or Aligned.sortedByCoord.out.patched.bam,
# or .bam (the GTEx bam directory has files with all 3 extensions, I'm not sure of the differences)

mkdir -p ${DIR_OUT}

echo "Running task $LSB_JOBINDEX"
echo "LIBID: $LIBID"

salmon quant \
  -i ${TXINDEX_salmon} --libType A \
  --validateMappings --seqBias \
  -p ${THREADS} --numBootstraps 50 \
  -1 ${DIR_FASTQ}/${LIBID}_R1.fastq \
  -2 ${DIR_FASTQ}/${LIBID}_R2.fastq \
  -o ${DIR_OUT}/${SAMPLE}