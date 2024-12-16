#!/usr/bin/env Rscript

# This script will perform marker protein analyses in a one-vs-all comparison for all manual_l2l3.

require(Seurat)
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

hc_pf_seuratobject_rds <- args[1]
findmarkers_list_rds <- args[2]

seuratObject <- readRDS(hc_pf_seuratobject_rds)

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(celltype = manual_l3,
                celltype = ifelse(celltype %in% c("Classical monocytes", "Macrophages C1Q+", "Macrophages SPP1+", "Macrophages VCAN+", "Macrophages VCAN+C1Q+"), "Mono-macs", as.character(celltype)))

Idents(seuratObject) <- "celltype"
seuratObject <- NormalizeData(seuratObject, assay = "CITE", normalization.method = "CLR")

markers_list <- FindAllMarkers(seuratObject, assay = "CITE", only.pos = T)
markers_list <- split(markers_list, markers_list$cluster)

saveRDS(markers_list, findmarkers_list_rds, compress = "gzip")

sessionInfo()
