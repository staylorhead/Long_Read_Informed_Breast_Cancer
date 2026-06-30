# summarize results on seadragon

library(data.table)
library(susieR)
library(ComplexUpset)
library(future.apply)
library(ComplexUpset)
library(dplyr)
library(DT)
library(tibble)

summarize_susie_results <- function(result, pip_thresh_8 = 0.8, pip_thresh_9 = 0.9) {
  if (is.null(result) || is.null(result$susie_fit)) return(NULL)
  
  susie_fit <- result$susie_fit
  pips <- susie_fit$pip
  snps <- colnames(susie_fit$alpha)
  
  names(pips) <- snps
  pips_sorted <- sort(pips, decreasing = TRUE)

  # Extract high PIPs
  pips_high_8 <- pips[pips > pip_thresh_8]
  pips_high_9 <- pips[pips > pip_thresh_9]
  
  # High-confidence credible sets from susie_fit$sets$cs
  cs <- susie_fit$sets$cs
  cs <- cs[!sapply(cs, is.null)]  # remove any NULL entries
  all_cs_snps <- unlist(cs, use.names = FALSE)
  all_cs_ids <- snps[all_cs_snps]
  all_cs_pips <- pips[all_cs_ids]
  
  summary <- tibble(
    transcript_id = result$transcript_id,
    num_cs = length(cs),
    cs_coverage = list(susie_fit$sets$coverage),
    top_pips = list(head(pips_sorted, 5)),
    high_pips_8 = list(pips_high_8),
    high_pips_9 = list(pips_high_9),
    cs_snps = list(all_cs_ids),
    cs_pips = list(all_cs_pips)
  )
  
  return(summary)
}

# summarize gene level
wd <- "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap"
file_list <- list.files(wd, full.names = TRUE)
#file_list <- file_list[-grep("aggregated",file_list)]

keep <- c(keep,grep("gene_gencodev45_k15",file_list))
keep <- c(keep,grep("gene_comb_k15",file_list))
keep <- c(keep,grep("gene_veigaTum_k200",file_list))
keep <- unique(keep)

file_list <- file_list[keep]
file_list <- file_list[grep("rmcov1",file_list)]
file_list <- file_list[grep("gene",file_list)]
file_list <- file_list[grep("_L",file_list)]
file_list <- file_list[grep("TCGA",file_list)]

#file_list <- file_list[grep("_L5",file_list)]
#keep <- c(keep,grep("tx_gencodev45_k15",file_list))
#keep <- c(keep,grep("tx_comb_k15",file_list))
#keep <- c(keep,grep("tx_veigaNorm_k50",file_list))
#keep <- unique(keep)

file_list <- file_list[file.info(file_list)$size > 0]

summary_results_list <- list()

for (i in seq_along(file_list)) {
  message("Processing file ", i, " of ", length(file_list))
  load(file_list[i])  # loads fullres

  set_name <- sub(".RData$", "", basename(file_list[i]))

  for (res in fullres) {
    tx_name <- as.character(res$transcript_id)
    if (!is.null(tx_name)) {
      summary_data <- summarize_susie_results(res)
      if (!is.null(summary_data)) {
        summary_data$set_name <- set_name
        summary_results_list[[length(summary_results_list) + 1]] <- summary_data
      }
    }
  }
}

summary_all <- bind_rows(summary_results_list)
# Combine all results efficiently
summary_results <- rbindlist(summary_results_list, use.names = TRUE, fill = TRUE)
save(summary_results,file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap/TCGA_aggregated_finemap_res_byL_gene.RData")




# summarize tx level
wd <- "/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap"
file_list <- list.files(wd, full.names = TRUE)
#file_list <- file_list[-grep("aggregated",file_list)]
keep <- {}
keep <- c(keep,grep("tx_gencodev45_k15",file_list))
keep <- c(keep,grep("tx_comb_k15",file_list))
keep <- c(keep,grep("tx_veigaTum_k15",file_list))
keep <- unique(keep)
file_list <- file_list[keep]
file_list <- file_list[grep("rmcov1",file_list)]
file_list <- file_list[grep("tx",file_list)]
file_list <- file_list[grep("_L",file_list)]
file_list <- file_list[grep("TCGA",file_list)]


#file_list <- file_list[grep("_L5",file_list)]
#keep <- c(keep,grep("tx_gencodev45_k15",file_list))
#keep <- c(keep,grep("tx_comb_k15",file_list))
#keep <- c(keep,grep("tx_veigaNorm_k50",file_list))
#keep <- unique(keep)

file_list <- file_list[file.info(file_list)$size > 0]

summary_results_list <- list()

for (i in seq_along(file_list)) {
  message("Processing file ", i, " of ", length(file_list))
  load(file_list[i])  # loads fullres

  set_name <- sub(".RData$", "", basename(file_list[i]))

  for (res in fullres) {
    tx_name <- as.character(res$transcript_id)
    if (!is.null(tx_name)) {
      summary_data <- summarize_susie_results(res)
      if (!is.null(summary_data)) {
        summary_data$set_name <- set_name
        summary_results_list[[length(summary_results_list) + 1]] <- summary_data
      }
    }
  }
}

summary_all <- bind_rows(summary_results_list)
# Combine all results efficiently
summary_results <- rbindlist(summary_results_list, use.names = TRUE, fill = TRUE)
save(summary_results,file="/rsrch5/home/epi/sthead/isoqtl_lr_breast/results/finemap/TCGA_aggregated_finemap_res_byL_tx.RData")


