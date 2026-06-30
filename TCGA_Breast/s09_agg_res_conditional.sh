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
num_perm=${12}

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools/TCGA_cond_${prefix}"

echo "Processing $dir"

cd $dir

# Count and concatenate normalized TUMOR files
count_tum_norm=$(ls conditional_gene_tumor_nperm${num_perm}_normalized_[0-9]*.txt | wc -l)
if [ "$count_tum_norm" -eq 22 ]; then
echo "  Found 22 TUM normalized files — combining..."
cat conditional_gene_tumor_nperm${num_perm}_normalized_[0-9]*.txt | gzip -c > gene_qtls_cond_nperm${num_perm}_tum_normalized_full.txt.gz
else
echo "  Skipping TUM normalized (found $count_tum_norm files)"
fi

# Count and concatenate normalized TUMOR files
count_tum_norm=$(ls conditional_tx_tumor_nperm${num_perm}_normalized_[0-9]*.txt | wc -l)
if [ "$count_tum_norm" -eq 22 ]; then
echo "  Found 22 TUM normalized files — combining..."
cat conditional_tx_tumor_nperm${num_perm}_normalized_[0-9]*.txt | gzip -c > tx_qtls_cond_nperm${num_perm}_tum_normalized_full.txt.gz
else
echo "  Skipping TUM normalized (found $count_tum_norm files)"
fi

# unload modules
module unload R


echo "***** job ends ***** "
date
