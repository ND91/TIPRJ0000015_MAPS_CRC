#!/usr/bin/env Rscript

# This script will perform marker gene expression analyses comparing the celltypes in TX.

require(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_rds <- args[2]
level <- args[3]

seuratObject <- readRDS(seurat_rds)

Idents(seuratObject) <- level

marker_genes <- FindAllMarkers(seuratObject, only.pos = F)
marker_genes <- split(marker_genes, f = marker_genes$cluster)

saveRDS(marker_genes, markers_rds, compress = "gzip")

sessionInfo()