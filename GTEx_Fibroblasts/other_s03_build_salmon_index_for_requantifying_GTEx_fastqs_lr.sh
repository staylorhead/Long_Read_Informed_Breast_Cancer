#BSUB -W 8:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q medium

# cd /rsrch5/home/epi/sthead/isoqtl_lr_breast

# bsub -J "salmon_idx" -W 8:00 -n 10 -M 80 -R rusage[mem=80] -q medium \
# -o /rsrch5/home/epi/sthead/isoqtl_lr_breast/logs/other/build_salmon_index_for_requant_GTEx_fibro_lr_%J.log \
# /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/GTEx_fibro/other_s03_build_salmon_index_for_requantifying_GTEx_fastqs_lr.sh

GENOME=/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta
GTF=/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/Cells_Cultured_fibroblasts_filtered/Cells_Cultured_fibroblasts_filtered_cleaned.gtf
DIR_OUT=/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v9/SQANTI3/orthogonal/Expression/short-read/salmon/GTEx_fibroblast/ESPRESSO_fibroblast_filtered

mkdir -p ${DIR_OUT}

source /rsrch5/home/epi/bhattacharya_lab/software/gffread/bin/activate

gffread ${GTF} \
  -g ${GENOME} \
  -w ${DIR_OUT}/GCA_000001405.15_GRCh38_no_alt_analysis_set_ESPRESSO_fibroblast_filtered.fasta
  
source /rsrch5/home/epi/bhattacharya_lab/software/gffread/bin/deactivate

source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/activate

sh /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/generateDecoyTranscriptome.sh \
 -j 10 \
 -b /risapps/rhel8/miniforge3/24.5.0-0/envs/bedtools-2.30.0/bin/bedtools \
 -a ${GTF} \
 -g ${GENOME} \
 -t ${DIR_OUT}/GCA_000001405.15_GRCh38_no_alt_analysis_set_ESPRESSO_fibroblast_filtered.fasta \
 -o ${DIR_OUT}
 
source /rsrch5/home/epi/bhattacharya_lab/software/mashmap/bin/deactivate

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

conda activate salmon-1.10.2

salmon index -t ${DIR_OUT}/gentrome.fa \
  -i ${DIR_OUT}/index \
  --decoys ${DIR_OUT}/decoys.txt \
  -k 31