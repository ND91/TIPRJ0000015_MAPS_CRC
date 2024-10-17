#!/usr/bin/env Rscript
# The goal of this script is to subset the kidney-derived macrophages from adults only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

library(Seurat)
library(dplyr)

seurat_rds <- args[1]
hc_kidney_macrophages_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(tissue_original %in% c("AdultKidney"),
                manual_l2 %in% c("Macrophages"))

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected$CellID)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 7983251)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:21)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:21, seed.use = 632632)

# Save data
saveRDS(seuratObject, hc_kidney_macrophages_seurat_rds, compress = "gzip")

sessionInfo()