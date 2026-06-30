library(data.table)
library(future)


# read in and aggregate data on seadragon
file_list <- list.files("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc")
file_list <- file_list[grep("GTEx",file_list)]
file_list <- file_list[-grep("fibro",file_list)]

#plan(multisession, workers = 6)  # adjust based on your system

# Define the function that processes one file
process_file <- function(file) {
  load(paste0("/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc/", file))
  set_name <- sub(".RData", "", file)
  
  # Get summary stats
  summary_df <- rbindlist(
    lapply(names(coloc_results), function(gene_id) {
      s <- coloc_results[[gene_id]]$summary
      as.data.table(as.list(s))[, gene := gene_id]
    }),
    fill = TRUE
  )
  summary_df[, set_name := set_name]
  
  # Get top SNPs
  top_snp_df <- rbindlist(
    lapply(names(coloc_results), function(gene_id) {
      res_df <- coloc_results[[gene_id]]$results
      if (!is.null(res_df) && "SNP.PP.H4" %in% colnames(res_df)) {
        top_row <- res_df[which.max(res_df$SNP.PP.H4), ]
        data.table(
          gene = gene_id,
          top_snp = top_row$snp,
          top_snp_pph4 = top_row$SNP.PP.H4
        )
      } else {
        data.table(gene = gene_id, top_snp = NA, top_snp_pph4 = NA)
      }
    }),
    fill = TRUE
  )
  
  merge(summary_df, top_snp_df, by = "gene", all.x = TRUE)
}

bad_files <- c()
res <- lapply(file_list, function(file) {
  tryCatch({
    process_file(file)
  }, error = function(e) {
    message("Error in file: ", file, " -> ", e$message)
    bad_files <<- c(bad_files, file)
    NULL
  })
})

res_all <- rbindlist(res, fill = TRUE)

save(res_all, file = "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/coloc/GTEx_aggregate_res.RData")