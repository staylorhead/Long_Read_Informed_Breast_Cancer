require(data.table)
require(isotwas)

trait <- c("HER2enriched","TN","lumA","lumBHER2neg","overallBC","lumB","PMD","DA","NDA","BCsurvival")


base_path="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged"

dirs <- list.dirs("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/isoTWAS/TCGA/output_overlap_ss_munged")
dirs <- dirs[grep("associations",dirs)]

for (quan in c("gencodev45", "veigaTum", "comb")) {
    print(quan)
    cor_dir <- dirs[grep(quan,dirs)]

    for (cancer in trait) {
    
    file_path <- file.path(cor_dir, "SignificantAssociations_isoTWAS_noFineMap.tsv")
    
    if (file.exists(file_path)) {
      sig_tx <- fread(file_path,fill=T)
      sig_tx$Quant <- quan

      out_file <- file.path(base_path, "isoTWAS_noFM_overlap_ss_munged.txt")
      fwrite(sig_tx, out_file, row.names = FALSE, quote = FALSE, append = TRUE,sep="\t")  # append if combining results
      cat("✓ Processed:", cancer, quan, "\n")
    } else {
      cat(" File not found:", file_path, "\n")
    }
  }
}
