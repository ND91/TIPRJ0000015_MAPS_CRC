#!/usr/bin/env Rscript

# This script will perform marker gene expression analyses comparing the different macrophages in TX.

require(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

Idents(seuratObject) <- "manual_l4"

l4_markers <- FindAllMarkers(seuratObject, only.pos = F)
l4_markers <- split(l4_markers, f = l4_markers$cluster)

saveRDS(l4_markers, markers_rds, compress = "gzip")

sessionInfo()