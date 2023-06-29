#!/usr/bin/env Rscript

# The goal of this script is to select the TX-derived cells from the CRC patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crc_tx_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

txcells <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX") %>%
  dplyr::pull(cellID)

seuratObject <- seuratObject[,seuratObject@meta.data$cellID %in% txcells]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 1432)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:68)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:68, seed.use = 631623)

# Save data
saveRDS(seuratObject, crc_tx_seurat_rds, compress = "gzip")

sessionInfo()