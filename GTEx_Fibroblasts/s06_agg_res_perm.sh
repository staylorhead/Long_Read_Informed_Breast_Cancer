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

dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools/GTEx_fibro_perm_${prefix}"

echo "Processing $dir"

cd $dir

# Count and concatenate normalized gene files
count_tumor_norm=$(ls gene_qtls_perm${num_perm}_normalized_[0-9]*.txt | wc -l)
if [ "$count_tumor_norm" -eq 22 ]; then
echo "  Found 22 normalized gene files — combining..."
cat gene_qtls_perm${num_perm}_normalized_*.txt | gzip -c > gene_qtls_perm${num_perm}_normalized_full.txt.gz
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/child_qtltools_runFDR_cis.R gene_qtls_perm${num_perm}_normalized_full.txt.gz 0.2 agg_gene_qtls_perm${num_perm}_normalized_full_FDR0.2.results
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/child_qtltools_runFDR_cis.R gene_qtls_perm${num_perm}_normalized_full.txt.gz 0.05 agg_gene_qtls_perm${num_perm}_normalized_full_FDR0.05.results
else
echo "  Skipping gene normalized (found $count_tumor_norm files)"
fi

# Count and concatenate normalized tx files
count_tumor_norm=$(ls tx_qtls_perm${num_perm}_normalized_[0-9]*.txt | wc -l)
if [ "$count_tumor_norm" -eq 22 ]; then
echo "  Found 22 normalized tx files — combining..."
cat tx_qtls_perm${num_perm}_normalized_*.txt | gzip -c > tx_qtls_perm${num_perm}_normalized_full.txt.gz
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/child_qtltools_runFDR_cis.R tx_qtls_perm${num_perm}_normalized_full.txt.gz 0.2 agg_tx_qtls_perm${num_perm}_normalized_full_FDR0.2.results
Rscript /rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/child_qtltools_runFDR_cis.R tx_qtls_perm${num_perm}_normalized_full.txt.gz 0.05 agg_tx_qtls_perm${num_perm}_normalized_full_FDR0.05.results
else
echo "  Skipping tx normalized (found $count_tumor_norm files)"
fi


# unload modules
module unload R


echo "***** job ends ***** "
date
