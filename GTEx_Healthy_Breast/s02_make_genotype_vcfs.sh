module load plink/2.00-alpha
module load tabix/0.2.6
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

# navigate to directory
cd /rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/geno/GTEx_MAF0.01_passQC

# create vcf for unrelated samples (none to remove here)

plink2 --pfile allchr \
--chr 1-22 \
--export vcf bgz id-paste=iid \
--out autosome_unrelated

# tabix

tabix -p vcf autosome_unrelated.vcf.gz

plink2 --vcf autosome_unrelated.vcf.gz \
--freq \
--out autosome_unrelated_vcf_MAF

# sort vcf file
bcftools sort autosome_unrelated.vcf.gz -Oz -o autosome_unrelated_sorted.vcf.gz

# add "chr" prefix to all rows in vcf
bcftools annotate \
  --rename-chrs /rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/bed_map_file \
  -Oz \
  -o autosome_unrelated_sorted_chr_labels.vcf.gz \
  autosome_unrelated_sorted.vcf.gz
  
# tabix new file
tabix -p vcf autosome_unrelated_sorted_chr_labels.vcf.gz