#!/usr/bin/env Rscript

# This script will perform the TAM classification analyses of PF and TX samples from paired donors by using UCell to predict TAM subclasses as defined by Ma et al. 2022 (DOI: https://doi.org/10.1016/j.it.2022.04.008).

BiocManager::install("UCell")

require(Seurat)
require(dplyr)
require(UCell)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

crcpmp_pf_tx_paired_macrophages_seuratobject_rds <- args[1]
tam_markers_xlsx <- args[2]
threads <- args[3]
crcpmp_pf_tx_paired_macrophages_tamannotation_seuratobject_rds <- args[4]

seuratObject <- readRDS(crcpmp_pf_tx_paired_macrophages_seuratobject_rds)

tam_markers <- readxl::read_excel(tam_markers_xlsx)
tam_markers_list <- lapply(split(tam_markers, f = tam_markers$TAM), function(tam_genes){
  tam_genes[(tam_genes$Gene %in% rownames(GetAssayData(seuratObject, assay = "RNA"))), ]$Gene
})

DefaultAssay(seuratObject) <- "RNA"
seuratObject_annotated <- AddModuleScore_UCell(seuratObject, 
                                               features = tam_markers_list, 
                                               name = NULL,
                                               ncores = threads)

saveRDS(seuratObject_annotated@meta.data, file = crcpmp_pf_tx_paired_macrophages_tamannotation_seuratobject_rds, compress = "gzip")

sessionInfo()