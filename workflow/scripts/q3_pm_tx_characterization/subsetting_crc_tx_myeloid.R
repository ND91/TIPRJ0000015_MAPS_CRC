#!/usr/bin/env Rscript

# The goal of this script is to select the TX-derived myeloid cells from the CRC patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds_path <- args[1]
crc_tx_myeloid_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

txmyeloidcells <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX",
                manual_l1 == "Myeloid") %>%
  dplyr::pull(cellID)

seuratObject <- seuratObject[,seuratObject@meta.data$cellID %in% txmyeloidcells]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 96781243)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:71)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:71, seed.use = 6132)

# Save data
saveRDS(seuratObject, crc_tx_myeloid_seurat_rds_path, compress = "gzip")

sessionInfo()