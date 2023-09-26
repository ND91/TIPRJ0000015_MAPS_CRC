#!/usr/bin/env Rscript

# The goal of this script is to select the PF macrophages from HC and CRC PM+ patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_crcpmp_pf_macrophages_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Group %in% c("CRC+", "HC"),
                Tissue %in% c("PF"),
                manual_l2 %in% c("Macrophages")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 4132)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:47)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:47, seed.use = 5132)

# Save data
saveRDS(seuratObject, hc_crcpmp_pf_macrophages_seurat_rds, compress = "gzip")

sessionInfo()