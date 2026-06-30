

module load plink/2.00-alpha
module load tabix/0.2.6
module load bcftools

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate bcftools-1.16

# navigate to directory with 1KG imputed variants
cd /rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/AffySNP6/VCF/normal/imputed

# perform QC of bcf file
mkdir -p MAF0.01_passQC

plink2 --bcf TCGA_BRCA_AffySNP6_phased_1kGP_imputed.bcf dosage=HDS \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.05 \
--make-king-table \
--make-pgen \
--out MAF0.01_passQC/allchr

# rename log file	
mv MAF0.01_passQC/allchr.log MAF0.01_passQC/allchr_variant_filter.log

# calculate variant- and sample-level missingness with filtered variants
plink2 --pfile MAF0.01_passQC/allchr \
--missing \
--out MAF0.01_passQC/allchr

# rename log file	
mv MAF0.01_passQC/allchr.log MAF0.01_passQC/allchr_missingness.log

# extract imputation quality scores from BCF file
bcftools view -h TCGA_BRCA_AffySNP6_phased_1kGP_imputed.bcf | grep "##INFO="

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/R2\t%INFO/ER2\t%INFO/AVG_CS\t%INFO/IMPUTED\t%INFO/TYPED\n' TCGA_BRCA_AffySNP6_phased_1kGP_imputed.bcf > imputation_quality.txt
gzip imputation_quality.txt

# filter bcf to include only those variants with imputation R2 > 0.8
bcftools filter -O z -o TCGA_BRCA_AffySNP6_phased_1kGP_imputed_impR2_08.bcf -i 'INFO/R2>0.8' TCGA_BRCA_AffySNP6_phased_1kGP_imputed.bcf

# perform QC of filtered bcf file
mkdir -p MAF0.01_passQC_impR2_08

plink2 --bcf TCGA_BRCA_AffySNP6_phased_1kGP_imputed_impR2_08.bcf dosage=HDS \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.05 \
--make-king-table \
--make-pgen \
--out MAF0.01_passQC_impR2_08/allchr
	
# rename log file	
mv MAF0.01_passQC_impR2_08/allchr.log MAF0.01_passQC_impR2_08/allchr_variant_filter.log

# calculate variant- and sample-level missingness with filtered variants
plink2 --pfile MAF0.01_passQC_impR2_08/allchr \
--missing \
--out MAF0.01_passQC_impR2_08/allchr

# rename log file	
mv MAF0.01_passQC_impR2_08/allchr.log MAF0.01_passQC_impR2_08/allchr_missingness.log


# perform QC of 0.3 filtered bcf file
mkdir -p MAF0.01_passQC_impR2_03

plink2 --vcf TCGA_BRCA_AffySNP6_phased_1kGP_imputed_impR2_03.bcf.gz dosage=HDS \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.05 \
--make-king-table \
--make-pgen \
--out MAF0.01_passQC_impR2_03/allchr
	
# rename log file	
mv MAF0.01_passQC_impR2_03/allchr.log MAF0.01_passQC_impR2_03/allchr_variant_filter.log

# calculate variant- and sample-level missingness with filtered variants
plink2 --pfile MAF0.01_passQC_impR2_03/allchr \
--missing \
--out MAF0.01_passQC_impR2_03/allchr

# rename log file	
mv MAF0.01_passQC_impR2_03/allchr.log MAF0.01_passQC_impR2_03/allchr_missingness.log