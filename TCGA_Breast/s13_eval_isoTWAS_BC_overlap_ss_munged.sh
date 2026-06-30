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

bin_index=${LSB_JOBINDEX}
isotwas_folder=$1
ld_folder=$2
gwas_file=$3
label=$4
trait=$5

module load R/4.3.1
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/TCGA/child_eval_isoTWAS_BC_overlap_ss_munged.R ${bin_index} ${isotwas_folder} ${ld_folder} ${gwas_file} ${label} ${trait}