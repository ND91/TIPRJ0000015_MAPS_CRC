#!/usr/bin/env Rscript

# This script will perform the TAM classification analyses by using UCell to predict TAM subclasses as defined by Ma et al. 2022 (DOI: https://doi.org/10.1016/j.it.2022.04.008).

BiocManager::install("UCell")

require(Seurat)
require(dplyr)
require(UCell)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

hc_pf_macrophages_seuratobject_rds <- args[1]
tam_markers_xlsx <- args[2]
threads <- args[3]
hc_pf_macrophages_tamannotation_responses_rds <- args[4]
hc_pf_macrophages_tamannotation_ranks_rds <- args[5]

seuratObject <- readRDS(hc_pf_macrophages_seuratobject_rds)

tam_markers <- readxl::read_excel(tam_markers_xlsx)
tam_markers_list <- lapply(split(tam_markers, f = tam_markers$TAM), function(tam_genes){
  tam_genes[(tam_genes$Gene %in% rownames(GetAssayData(seuratObject, assay = "RNA"))), ]$Gene
})

DefaultAssay(seuratObject) <- "RNA"
seuratObject_annotated <- AddModuleScore_UCell(seuratObject, 
                                               features = tam_markers_list, 
                                               name = NULL,
                                               ncores = threads)

cell_tam_markers_ranks <- StoreRankings_UCell(GetAssayData(seuratObject, assay = "RNA"))

saveRDS(seuratObject_annotated@meta.data, file = hc_pf_macrophages_tamannotation_responses_rds, compress = "gzip")
saveRDS(cell_tam_markers_ranks, file = hc_pf_macrophages_tamannotation_ranks_rds, compress = "gzip")

sessionInfo()