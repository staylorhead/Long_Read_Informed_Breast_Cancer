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

exp_file=$1
gen_file=$2
cov_file=$3
label=$4
gtf_file=$5
chr=$6
index=${LSB_JOBINDEX}

module load R/4.3.1

Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/GTEx_fibro/child_fit_isoTWAS_model_overlap_ss_munged.R ${exp_file} ${gen_file} ${cov_file} ${label} ${gtf_file} ${chr} ${index}
