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
gen_pcs_file=${11}

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

# script below will generate files for both nat and tumor and gene and tx if they don't already exist
cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/GTEx_cov_qtltools_gene_tumor_${prefix}.txt"

if [[ ! -f "${cov_file}" ]]; then
    echo "Generating HCP covariate files..."
    Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/GTEx/child_HCP_by_param.R ${assembly} ${do_ods} ${lstpm} ${do_vst} ${do_tmm} ${do_bc} ${k} ${method} ${gen_pcs_file} ${num_gen_pcs}
else
    echo "Covariate files already exist..."
fi

# load other modules
module load qtltools
module load tabix
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_gene_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}.bed"

if [[ ! -f "${bed_file}.gz" ]]; then
    echo "bgzipping and tabixing gene bed file..."
    bgzip ${bed_file} && tabix -p bed ${bed_file}.gz
else
    echo "bed file already bgzipped and tabixed..."
fi

bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/GTEx_tx_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}.bed"

if [[ ! -f "${bed_file}.gz" ]]; then
    echo "bgzipping and tabixing tx bed file..."
    bgzip ${bed_file} && tabix -p bed ${bed_file}.gz
else
    echo "bed file already bgzipped and tabixed..."
fi
