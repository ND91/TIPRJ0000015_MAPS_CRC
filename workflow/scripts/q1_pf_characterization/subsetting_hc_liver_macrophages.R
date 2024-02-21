#!/usr/bin/env Rscript

# The goal of this script is to select the macrophages from HC Colon patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_pf_liver_colon_macrophages_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Tissue %in% c("Liver"),
                manual_l2 %in% c("Macrophages")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 743)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:46)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:46, seed.use = 74327)

# Save data
saveRDS(seuratObject, hc_pf_liver_colon_macrophages_seurat_rds, compress = "gzip")

sessionInfo()