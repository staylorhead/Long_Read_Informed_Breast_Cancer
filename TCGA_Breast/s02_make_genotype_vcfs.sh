module load plink/2.00-alpha
module load tabix/0.2.6
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

# navigate to directory
cd /rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC

# PCA

plink2 --pfile allchr \
--keep /rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt \
--indep-pairwise 50 10 0.1 \
--autosome \
--out indep_pairwise_snps_in_unrelated
# 2428455/2904977 variants removed

plink2 --pfile allchr \
--keep /rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt \
--extract indep_pairwise_snps_in_unrelated.prune.in \
--write-samples \
--pca 50 \
--out autosome_unrelated

# make new genotype files with participant ID
awk 'NR>1 {print $1, substr($1, 1, 12)}' autosome_unrelated.id > autosome_unrelated_participant.id

# create pgen for unrelated samples with full ID
plink2 --pfile allchr \
--keep /rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt \
--autosome \
--make-pgen \
--out autosome_unrelated_full_ids
# 2904977 out of 2908600 variants loaded from allchr.pvar

# create vcf for unrelated samples with short (participant) ID
plink2 --pfile autosome_unrelated_full_ids \
--update-ids autosome_unrelated_participant.id \
--export vcf bgz vcf-dosage=DS-force id-paste=iid \
--out autosome_unrelated_participant_ids_GT_and_DS

plink2 --pfile autosome_unrelated_full_ids \
--update-ids autosome_unrelated_participant.id \
--export vcf bgz id-paste=iid \
--out autosome_unrelated_participant_ids_GT

# tabix

tabix -p vcf autosome_unrelated_participant_ids_GT.vcf.gz
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz

plink2 --vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz dosage=DS \
--freq \
--out autosome_unrelated_participant_ids_vcf_DS_MAF

plink2 --vcf autosome_unrelated_participant_ids_GT.vcf.gz \
--freq \
--out autosome_unrelated_participant_ids_vcf_GT_MAF

module load bcftools
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

bcftools sort autosome_unrelated_participant_ids_GT.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_sorted.vcf.gz

bcftools sort autosome_unrelated_participant_ids_GT_and_DS.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz

bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_sorted.vcf.gz
  
bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz
  
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz

tabix -p vcf autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz








# Now, repeat for imputation filtered data

cd /rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_08

# create pgen for unrelated samples with full ID
plink2 --pfile allchr \
--keep /rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt \
--autosome \
--make-pgen \
--out autosome_unrelated_full_ids
# 1122511 out of 1122530 variants loaded from allchr.pvar

# create vcf for unrelated samples with short (participant) ID
plink2 --pfile autosome_unrelated_full_ids \
--update-ids ../MAF0.01_passQC/autosome_unrelated_participant.id \
--export vcf bgz vcf-dosage=DS-force id-paste=iid \
--out autosome_unrelated_participant_ids_GT_and_DS

plink2 --pfile autosome_unrelated_full_ids \
--update-ids ../MAF0.01_passQC/autosome_unrelated_participant.id \
--export vcf bgz id-paste=iid \
--out autosome_unrelated_participant_ids_GT

# tabix

tabix -p vcf autosome_unrelated_participant_ids_GT.vcf.gz
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz

plink2 --vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz dosage=DS \
--freq \
--out autosome_unrelated_participant_ids_vcf_DS_MAF

plink2 --vcf autosome_unrelated_participant_ids_GT.vcf.gz \
--freq \
--out autosome_unrelated_participant_ids_vcf_GT_MAF

module load bcftools
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

bcftools sort autosome_unrelated_participant_ids_GT.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_sorted.vcf.gz

bcftools sort autosome_unrelated_participant_ids_GT_and_DS.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz

bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_sorted.vcf.gz
  
bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz
  
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz

tabix -p vcf autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz






# Now, repeat for imputation 0.3 filtered data

cd /rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed/MAF0.01_passQC_impR2_03

# create pgen for unrelated samples with full ID
plink2 --pfile allchr \
--keep /rsrch5/home/epi/sthead/isoqtl_TCGA/files_for_analysis/geno_sample_ids_unrelated_for_pca.txt \
--autosome \
--make-pgen \
--out autosome_unrelated_full_ids
# 2679319 out of 2681601 variants loaded from allchr.pvar.

# create vcf for unrelated samples with short (participant) ID
plink2 --pfile autosome_unrelated_full_ids \
--update-ids ../MAF0.01_passQC/autosome_unrelated_participant.id \
--export vcf bgz vcf-dosage=DS-force id-paste=iid \
--out autosome_unrelated_participant_ids_GT_and_DS

plink2 --pfile autosome_unrelated_full_ids \
--update-ids ../MAF0.01_passQC/autosome_unrelated_participant.id \
--export vcf bgz id-paste=iid \
--out autosome_unrelated_participant_ids_GT

# tabix

tabix -p vcf autosome_unrelated_participant_ids_GT.vcf.gz
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz

plink2 --vcf autosome_unrelated_participant_ids_GT_and_DS.vcf.gz dosage=DS \
--freq \
--out autosome_unrelated_participant_ids_vcf_DS_MAF

plink2 --vcf autosome_unrelated_participant_ids_GT.vcf.gz \
--freq \
--out autosome_unrelated_participant_ids_vcf_GT_MAF

module load bcftools
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

bcftools sort autosome_unrelated_participant_ids_GT.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_sorted.vcf.gz

bcftools sort autosome_unrelated_participant_ids_GT_and_DS.vcf.gz -Oz -o autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz

bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_sorted.vcf.gz
  
bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/TCGA/bed/bed_map_file \
  -Oz \
  -o autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz \
  autosome_unrelated_participant_ids_GT_and_DS_sorted.vcf.gz
  
tabix -p vcf autosome_unrelated_participant_ids_GT_and_DS_sorted_chr_labels.vcf.gz

tabix -p vcf autosome_unrelated_participant_ids_GT_sorted_chr_labels.vcf.gz
