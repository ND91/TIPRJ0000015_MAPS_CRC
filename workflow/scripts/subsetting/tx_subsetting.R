#!/usr/bin/env Rscript
# The goal of this script is to select the TX-derived cells only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
tx_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject <- seuratObject[,seuratObject@meta.data$Tissue == "TX"]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 12394870)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:63)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunTSNE(seuratObject, dims = 1:63, seed.use = 17298346)

# Save data
saveRDS(seuratObject, tx_seurat_rds_path, compress = "gzip")

sessionInfo()