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

# assembly="comb"
# do_ods=1
# lstpm=0
# do_vst=1
# do_tmm=1
# do_bc=0
# k=15
# method="HC"
# num_gen_pcs=5
# geno_file="/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03/autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz"

# load other modules
module load qtltools
module unload R
module load R/4.3.1

#### gene-level ######################################################

echo "BEGINNING GENE-LEVEL ANALYSIS..."

prefix="${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_k${k}_${method}_genpc${num_gen_pcs}"

cd /rsrch5/home/epi/sthead/isoqtl_lr_breast

out_dir="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/qtltools/TCGA_nom_${prefix}"
bed_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed/TCGA_gene_exp_${assembly}_ODS${do_ods}_TMM${do_tmm}_VST${do_vst}_lstpm${lstpm}_doBC${do_bc}_tumor.bed.gz"
cov_file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/covariates/TCGA_cov_qtltools_gene_tumor_${prefix}.txt"

mkdir -p ${out_dir}
cd ${out_dir}

echo "Running QTLtools to get eQTL SS..."

QTLtools cis --vcf ${geno_file} --bed ${bed_file} --nominal 1 --chunk ${i} 100 --cov ${cov_file} --out nom_gene_chunk${i}.txt --normal --std-err

eqtlFile=${out_dir}/nom_gene_chunk${i}.txt

script="/rsrch5/home/epi/sthead/isoqtl_lr_breast/scripts/TCGA/child_run_coloc_BC.R"

echo "Running coloc..."

if [ -e "$eqtlFile" ]; then

  for pheno in overallBC lumA lumB lumBHER2neg HER2enriched TN; do

    echo $pheno

    level="gene"

    # set GWAS file path conditionally
    if [ "$pheno" == "BCsurvival" ]; then
      gwasFile=/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_Survival_BCAC2021/munged_GRCh38_survival_allpatients_ambig_strand_removed.tsv.gz
    else
      gwasFile=/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/Breast_Cancer_BCAC2020/munged_GRCh38_${pheno}_ambig_strand_removed.tsv.gz
    fi

    outfile=/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc/TCGA_${assembly}_${pheno}_${level}_results_chunk${i}.RData

    if [ ! -e "$outfile" ]; then
      echo "Running coloc for ${pheno}, chunk ${i}"
      Rscript "$script" "$gwasFile" "$eqtlFile" "$pheno" "$i" "$level" "$assembly"
    else
      echo "Skipping ${pheno}, chunk ${i} (output exists)"
    fi

  done

  # remove nominal qtltools results
  rm "$eqtlFile"

fi



