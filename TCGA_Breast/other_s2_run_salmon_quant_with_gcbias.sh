#!/bin/sh
#BSUB -J salmon_quant[1-2]
#BSUB -o /rsrch5/home/epi/sthead/isoqtl_TCGA/logs/quant_gcbias_%J_%I.out
#BSUB -e /rsrch5/home/epi/sthead/isoqtl_TCGA/logs/quant_gcbias_%J_%I.err
#BSUB -W 3:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q short

DIR_FASTQ=/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/fastq
DIR_OUT=/rsrch5/home/epi/bhattacharya_lab/data/Veiga/STB/SQANTI3/orthogonal/short-read/salmon/combined_gencodev45_and_ESPRESSO_tumor_filtered
TXINDEX_salmon=${DIR_OUT}/index
THREADS=10

INPUT_FILE="/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/fastq_hash_codes_to_TCGA_IDs.txt"
LIBID=$(sed -n "${LSB_JOBINDEX}p" "$INPUT_FILE" | cut -f1)

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

mkdir -p ${DIR_OUT}/run

conda activate salmon-1.10.2 

echo "Running task $LSB_JOBINDEX"
echo "LIBID: $LIBID"

salmon quant \
  -i ${TXINDEX_salmon} --libType A \
  --validateMappings --seqBias --gcBias \
  -p ${THREADS} --numBootstraps 50 \
  -1 ${DIR_FASTQ}/${LIBID}.rna_seq.transcriptome.gdc_realn.bam_R1_001.fastq.gz \
  -2 ${DIR_FASTQ}/${LIBID}.rna_seq.transcriptome.gdc_realn.bam_R2_001.fastq.gz \
  -o ${DIR_OUT}/run_with_gcbias_corr/${LIBID}
  