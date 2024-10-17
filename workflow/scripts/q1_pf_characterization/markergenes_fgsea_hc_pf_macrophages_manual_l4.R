#!/usr/bin/env Rscript

# This script will perform the marker gene pathway analyses analyses comparing the macrophage subsets at manual_l4 level with one another.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)
require(fgsea)
require(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

markergenes_list_rds <- args[1]
fgsea_list_rds <- args[2]
fgsea_list_xlsx <- args[3]

markergenes_list <- readRDS(markergenes_list_rds)

gs_hs <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

gs_hs_list <- split(x = gs_hs$gene_symbol, f = gs_hs$gs_name)

fgsea_list <- lapply(markergenes_list, function(macrophage){
  degs_df <- macrophage %>%
    data.frame() %>%
    dplyr::filter(!is.na(gene))
  
  es <- degs_df$avg_log2FC
  names(es) <- degs_df$gene
  
  fgsea_obj <- fgsea(pathways = gs_hs_list, 
                     stats    = es,
                     minSize  = 15,
                     maxSize  = 500)
  
  fgsea_obj <- fgsea_obj[order(fgsea_obj$pval),]
  
  return(fgsea_obj)
})

saveRDS(fgsea_list, fgsea_list_rds)

write.xlsx(lapply(fgsea_list, function(fgsea_obj){
  data.frame(fgsea_obj[,1:7])
}), file = fgsea_list_xlsx)


sessionInfo()
