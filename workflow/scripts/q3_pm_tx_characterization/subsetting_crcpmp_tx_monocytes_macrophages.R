#!/usr/bin/env Rscript

# The goal of this script is to select the TX-derived macrophages from the CRC patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crc_tx_monocytes_macrophages_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

txmonomac <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX",
                manual_l2 %in% c("Monocytes", "Macrophages")) %>%
  dplyr::pull(cellID)

seuratObject <- seuratObject[,seuratObject@meta.data$cellID %in% txmonomac]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 12363216)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:50)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:50, seed.use = 2135325)

# seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:46)
# seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
# seuratObject <- RunUMAP(seuratObject, dims = 1:46, seed.use = 62136)

# Save data
saveRDS(seuratObject, crc_tx_monocytes_macrophages_seurat_rds, compress = "gzip")

sessionInfo()