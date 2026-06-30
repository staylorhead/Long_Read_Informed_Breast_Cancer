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
level=$2
k=$3
finemap_dir=$4
L=$5
chr=${LSB_JOBINDEX}

module load R/4.3.1

Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/GTEx_fibro/child_finemap_ld_gwas.R "${assembly}" "${level}" "${k}" "${finemap_dir}" "${chr}" "${L}"

module unload R 