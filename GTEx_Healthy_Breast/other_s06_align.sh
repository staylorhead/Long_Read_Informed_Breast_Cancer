#BSUB -J star_align[1-440]%25
#BSUB -W 24:00
#BSUB -o /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/align_%J_%I.out
#BSUB -e /rsrch5/home/epi/sthead/isoqtl_GTEx/logs/align_%J_%I.out
#BSUB -q medium 
#BSUB -n 4
#BSUB -M 50
#BSUB -R rusage[mem=50]

### DEFINE GLOBAL VARIABLES

DIR_BAM=/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam
SAMPLES_FILE="/rsrch5/home/epi/sthead/isoqtl_lr_breast/files_for_analysis/GTEx_breast_mammary_bam_files_with_sample_attrib.txt"
ROWID=$((LSB_JOBINDEX + 1))
LIBID=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f1)
SAMPLE=$(sed -n "${ROWID}p" "$SAMPLES_FILE" | cut -f2)
DIR_FASTQ="/rsrch5/scratch/epi/sthead/GTEx_v8/fastq"
fasta_file="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_file="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v45/gencode.v45.annotation.nochr.gtf"
star_database_v45_grch38="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/star_database_v45_grch38"
out_folder_alignment="/rsrch5/scratch/epi/sthead/GTEx_v8/bam"
file1=${DIR_FASTQ}/${LIBID}_R1.fastq
file2=${DIR_FASTQ}/${LIBID}_R2.fastq
out_folder_trim=/rsrch5/scratch/epi/sthead/GTEx_v8/trimmed
file1_trimmed=${out_folder_trim}/${SAMPLE}_R1_trimmed.fastq.gz
file2_trimmed=${out_folder_trim}/${SAMPLE}_R2_trimmed.fastq.gz
out_folder_fastqc=/rsrch5/scratch/epi/sthead/GTEx_v8/fastqc

### RUN CODE

# fastqc

module load fastqc
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate fastqc-0.11.9

fastqc "$file1" "$file2" -o "$out_folder_fastqc" -d "/rsrch5/scratch/epi/sthead/GTEx_b8/tmp_fastqc"

conda deactivate

# trimmomatic

module add trimmomatic
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate trimmomatic-0.39

trimmomatic PE -threads 1 "$file1" "$file2" \
${out_folder_trim}/${SAMPLE}_R1_trimmed.fastq \
${out_folder_trim}/${SAMPLE}_R1_unpaired.fastq \
${out_folder_trim}/${SAMPLE}_R2_trimmed.fastq \
${out_folder_trim}/${SAMPLE}_R2_unpaired.fastq \
LEADING:20 TRAILING:20 MINLEN:50

conda deactivate

gzip ${out_folder_trim}/${SAMPLE}_R1_trimmed.fastq
gzip ${out_folder_trim}/${SAMPLE}_R2_trimmed.fastq

# star alignment

module load star

if [ -d "$star_database_v45_grch38" ]; then
    echo "Folder exists"
else
STAR --runMode genomeGenerate \
--genomeFastaFiles "$fasta_file" \
--sjdbGTFfile "$gtf_file" \
--genomeDir "$star_database_v45_grch38" \
--genomeSAindexNbases 12
fi
STAR --genomeDir "$star_database_v45_grch38" \
--readFilesIn "$file1_trimmed" "$file2_trimmed"  \
--outFileNamePrefix "$out_folder_alignment"/${SAMPLE}.sortedbyCoord_ \
--outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNmax 5 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--readFilesCommand zcat 

module unload star

### CLEANUP 

rm ${out_folder_trim}/${SAMPLE}_R1_*.fastq*
rm ${out_folder_trim}/${SAMPLE}_R2_*.fastq*
