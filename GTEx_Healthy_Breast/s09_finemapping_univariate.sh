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
level=$2
k=$3
nperm=$4
maf_filt=$5
rm_cov=$6
nsets=$7
L=$8

echo "Running fine mapping..."
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/GTEx/child_finemapping_univariate.R ${level} ${assembly} ${k} ${nperm} ${maf_filt} ${rm_cov} ${nsets} ${LSB_JOBINDEX} 1 ${L}

# unload modules
module unload R

echo "***** job ends ***** "
date
