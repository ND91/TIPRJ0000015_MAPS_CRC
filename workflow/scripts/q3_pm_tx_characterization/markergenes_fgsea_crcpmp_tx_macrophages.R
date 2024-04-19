#!/usr/bin/env Rscript

# This script will perform marker gene set enrichment analyses on the marker genes when comparing the different macrophages in TX.

require(fgsea)
require(openxlsx)
require(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

markers_list_rds <- args[1]
fgsea_markers_list_rds <- args[2]
fgsea_markers_list_xlsx <- args[3]

gs_hs <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
gs_hs_list <- split(x = gs_hs$gene_symbol, f = gs_hs$gs_name)

marker_genes_list <- readRDS(markers_list_rds)

fgsea_marker_list <- lapply(marker_genes_list, function(marker_genes){
  marker_genes <- marker_genes %>%
    data.frame() %>%
    dplyr::filter(!is.na(gene))
  
  stats <- marker_genes$avg_log2FC
  names(stats) <- marker_genes$gene
  
  fgsea_obj <- fgsea(pathways = gs_hs_list, 
                     stats    = stats,
                     minSize  = 15,
                     maxSize  = 500)
  
  fgsea_obj <- fgsea_obj[order(fgsea_obj$pval),]
  
  return(fgsea_obj)
})

saveRDS(fgsea_marker_list, fgsea_markers_list_rds, compress = "gzip")

write.xlsx(lapply(fgsea_marker_list, function(fgsea_obj){
  data.frame(fgsea_obj[,1:7])
}), file = fgsea_markers_list_xlsx)

sessionInfo()