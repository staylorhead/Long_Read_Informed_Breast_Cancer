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

ref=$1
target=$2
pheno=$3

echo "Running gwas formatting..."
Rscript /rsrch5/home/epi/bhattacharya_lab/users/sthead/isoqtl_TCGA/scripts/liftover_gwas_ss_BC.R ${ref} ${target} ${pheno}

# unload modules
module unload R

echo "***** job ends ***** "
date
