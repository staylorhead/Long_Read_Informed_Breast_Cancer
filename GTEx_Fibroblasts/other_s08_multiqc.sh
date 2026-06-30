

# Load or activate MultiQC if needed
module load multiqc  # or your preferred environment

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate multiqc-1.25.2

# Run MultiQC
multiqc /rsrch5/scratch/epi/sthead/GTEx_v8/picard/ \
  -o /rsrch5/scratch/epi/sthead/GTEx_v8/picard_multiqc_with_fibro/