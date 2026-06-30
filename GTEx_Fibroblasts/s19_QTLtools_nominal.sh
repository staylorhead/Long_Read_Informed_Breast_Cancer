#!/bin/bash

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "User: $LSB_USER"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

assembly=$1
do_ods=$2
lstpm=$3
do_vst=$4
do_tmm=$5
do_bc=$6
k=$7
method=$8
num_gen_pcs=$9
geno_file=${10}
i=${LSB_JOBINDEX}

# assembly="gencodev45"
# do_ods=1
# lstpm=0
# do_vst=1
# do_tmm=1
# do_bc=0
# k=15
# method="HC"
# num_gen_pcs=5
# geno_file="/rsrch5/home/epi/bhattacharya_lab/users/sthead/isoqtl_GTEx/files_for_analysis/geno/MAF0.01_passQC/autosome_unrelated_sorted_chr_labels.vcf.gz"

# load other modules
module load qtltools
module unload R
module load R/4.3.1

#### gene-level ######################################################

echo "BEGINNING GENE-LEVEL ANALYSIS..."

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

cd /rsrch5/home/epi/sthead/isoqtl_lr_breast

out_dir="/rsrch5/home/epi/bhattacharya_lab/users/sthead/isoqtl_lr_breast/results/qtltools/GTEx_fibro_nom_${prefix}"
bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_fibroblast_gene_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}.bed.gz"
cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_fibro_cov_qtltools_gene_${prefix}.txt"

mkdir -p ${out_dir}
cd ${out_dir}

echo "Running QTLtools to get eQTL SS..."

QTLtools cis --vcf ${geno_file} --bed ${bed_file} --nominal 1 --chunk ${i} 100 --cov ${cov_file} --out nom_gene_chunk${i}.txt --normal --std-err

#### tx-level ######################################################

echo "BEGINNING TX-LEVEL ANALYSIS..."

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

cd /rsrch5/home/epi/sthead/isoqtl_lr_breast

out_dir="/rsrch5/home/epi/bhattacharya_lab/users/sthead/isoqtl_lr_breast/results/qtltools/GTEx_fibro_nom_${prefix}"
bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_fibroblast_tx_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}.bed.gz"
cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_fibro_cov_qtltools_tx_${prefix}.txt"

mkdir -p ${out_dir}
cd ${out_dir}

echo "Running QTLtools to get eQTL SS..."

QTLtools cis --vcf ${geno_file} --bed ${bed_file} --nominal 1 --chunk ${i} 100 --cov ${cov_file} --out nom_tx_chunk${i}.txt --normal --std-err



