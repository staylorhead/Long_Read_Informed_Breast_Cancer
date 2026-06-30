module load plink

cd /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8

OUT_DIR="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis"

mkdir -p ${OUT_DIR}/geno/GTEx_MAF0.01_passQC

plink2 --bfile /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.05 \
--make-king-table \
--make-pgen \
--out ${OUT_DIR}/geno/GTEx_MAF0.01_passQC/allchr

# - Remove variants with a minor allele frequency (MAF) $<1\%$ (1625 variants removed)
# - Remove variants with $>5\%$ missing genotype data (0 variants removed)
# - Remove variants failing Hardy-Weinberg equilibrium (HWE) at $p <$ 1e-6 (45728 variants removed)

# **1131422** variants remained after these filters. I proceed to compute relatedness statistics (KING relatedness coefficient) among these samples and output a KING kinship table. 

# calculate missingness:

module load plink/2.00-alpha
module load tabix/0.2.6
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

cd /rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/geno/GTEx_MAF0.01_passQC

# rename log file	
mv allchr.log allchr_variant_filter.log

# calculate variant- and sample-level missingness with filtered variants
plink2 --pfile allchr \
--missing \
--out allchr

# rename log file	
mv allchr.log allchr_missingness.log