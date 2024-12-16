#!/usr/bin/env Rscript

# This script will perform the marker gene pathway analyses analyses comparing the macrophage subsets at manual_l4 level with one another.

require(Seurat)
require(dplyr)
require(DESeq2)
require(openxlsx)
require(fgsea)
require(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

markergenes_list_rds <- args[1]
fgsea_kegg_list_rds <- args[2]
fgsea_kegg_list_xlsx <- args[3]
fgsea_gobp_list_rds <- args[4]
fgsea_gobp_list_xlsx <- args[5]

markergenes_list <- readRDS(markergenes_list_rds)

gs_hs_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
gs_hs_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

gs_hs_kegg_list <- split(x = gs_hs_kegg$gene_symbol, f = gs_hs_kegg$gs_name)
gs_hs_gobp_list <- split(x = gs_hs_gobp$gene_symbol, f = gs_hs_gobp$gs_name)

fgsea_markergenes <- function(macrophage, geneset_list){
  degs_df <- macrophage %>%
    data.frame() %>%
    dplyr::filter(!is.na(gene))
  
  es <- degs_df$avg_log2FC
  names(es) <- degs_df$gene
  
  fgsea_obj <- fgsea(pathways = geneset_list, 
                     stats = es,
                     minSize = 15,
                     maxSize = 500)
  
  fgsea_obj <- fgsea_obj[order(fgsea_obj$pval),]
  
  return(fgsea_obj)
}

fgsea_kegg_list <- lapply(markergenes_list, FUN = fgsea_markergenes, geneset_list = gs_hs_kegg_list)
fgsea_gobp_list <- lapply(markergenes_list, FUN = fgsea_markergenes, geneset_list = gs_hs_gobp_list)

saveRDS(fgsea_kegg_list, fgsea_kegg_list_rds)
saveRDS(fgsea_gobp_list, fgsea_gobp_list_rds)

write.xlsx(lapply(fgsea_kegg_list, function(fgsea_obj){
  data.frame(fgsea_obj[,1:7])
}), file = fgsea_kegg_list_xlsx)

write.xlsx(lapply(fgsea_gobp_list, function(fgsea_obj){
  data.frame(fgsea_obj[,1:7])
}), file = fgsea_gobp_list_xlsx)

sessionInfo()
