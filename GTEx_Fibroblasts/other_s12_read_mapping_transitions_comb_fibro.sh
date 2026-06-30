#!/bin/bash

SAMPLE=$LSB_JOBINDEX
GENE=$1

module load R/4.3.1

Rscript /rsrch5/scratch/epi/sthead/sankey/scripts/other_s12_read_mapping_transitions_comb_fibro.R "$SAMPLE" "$GENE" 

module unload R