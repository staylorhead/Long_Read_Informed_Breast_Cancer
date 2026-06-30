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

# load R
module load R/4.3.1

assembly=$1
do_ods=$2
lstpm=$3
do_vst=$4
do_tmm=$5
do_bc=$6
k=$7
method=$8
geno_file=$9
num_gen_pcs=${10}
chunk=${LSB_JOBINDEX}
gen_pcs_file=${11}
num_perm=${12}

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

# load other modules
module load qtltools
module load tabix
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

cd /rsrch5/home/epi/sthead/isoqtl_lr_brest

out_dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools/GTEx_fibro_perm_${prefix}"
bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_fibroblast_tx_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}.bed"
cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_fibro_cov_qtltools_tx_${prefix}.txt"

mkdir -p ${out_dir}
cd ${out_dir}

echo "Running normalization and QTLtool permutation for tx-level..."
QTLtools cis --vcf ${geno_file} --bed ${bed_file}.gz --permute ${num_perm} --region chr${chunk} --cov ${cov_file} --out tx_qtls_perm${num_perm}_normalized_${chunk}.txt --normal --seed 1219

# unload modules
module unload R
module unload tabix
module unload qtltools
module unload bcftools

echo "***** job ends ***** "
date
