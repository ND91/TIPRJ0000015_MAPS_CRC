#!/usr/bin/env Rscript

# The goal of this script is to select the PBMC-, PF-, and TX-derived myeloid cells from the CRC patients that provided all three samples.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds_path <- args[1]
crc_pbmc_pf_tx_myeloid_paired_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject <- seuratObject[,seuratObject@meta.data$manual_l1 == "Myeloid"]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 965321)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:25)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:25, seed.use = 61216)

# Save data
saveRDS(seuratObject, crc_pbmc_pf_tx_myeloid_paired_seurat_rds_path, compress = "gzip")

sessionInfo()