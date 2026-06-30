# Long_Read_informed_Breast_Cancer_Paper

This repository contains scripts for the analyses described in:

> **Improving isoform-level eQTL and integrative genetic analyses of breast cancer risk with long-read RNA transcript assemblies**  
> Head et al. (2026) [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.22.713514v2)

## Overview

We developed a framework for tissue-informed isoform-level eQTL mapping and TWAS+colocalization of breast cancer risk and leverages publicly available long-read RNA-seq assemblies to define tissue-specific transcriptomes. Gene- and isoform-level expression was quantified in three datasets using three annotation strategies each: standard GENCODE v45, a tissue-specific long-read-derived assembly, and a combined annotation merging both.

The pipeline proceeds from raw genotype and RNA-seq data through expression quantification, eQTL mapping, fine-mapping, colocalization with breast cancer GWAS, and isoTWAS association testing.

## Datasets and Pipelines

Scripts are organized by dataset:

### GTEx v8 Breast Mammary Tissue
**Scripts:** `GTEx_Healthy_Breast` | [`README`](GTEx_Healthy_Breast/README.md)

- **Samples:** Unrelated GTEx v8 breast mammary tissue donors
- **Long-read annotation:** Veiga et al. ESPRESSO normal assembly (`veigaNorm`)

### GTEx v8 Cultured Fibroblasts
**Scripts:** `GTEx_Fibroblasts` | [`README`](GTEx_Fibroblasts/README.md)

- **Samples:** Unrelated GTEx v8 cultured fibroblast donors
- **Long-read annotation:** GTEx v9 SQANTI3 ESPRESSO fibroblast assembly (`ESPRESSO`)

### TCGA BRCA Breast Tumor
**Scripts:** `TCGA_Breast_Tumor` | [`README`](TCGA_Breast_Tumor/README.md)

- **Samples:** Primary solid tumor samples from TCGA BRCA; genotypes from imputed Affymetrix SNP6 arrays
- **Long-read annotation:** Veiga et al. ESPRESSO tumor assembly (`veigaTum`)
