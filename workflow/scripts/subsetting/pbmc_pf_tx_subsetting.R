#!/usr/bin/env Rscript
# The goal of this script is to perform clustering and dimension reduction on all samples.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
pbmc_pf_tx_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 7693412)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:76)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunTSNE(seuratObject, dims = 1:76, seed.use = 9871234)

# Save data
saveRDS(seuratObject, pbmc_pf_tx_seurat_rds_path, compress = "gzip")

sessionInfo()