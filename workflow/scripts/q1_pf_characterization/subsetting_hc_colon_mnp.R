#!/usr/bin/env Rscript

# The goal of this script is to select the mnp from HC colon patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_colon_mnp_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Tissue %in% c("Colon"),
                manual_l2 %in% c("Monocytes", "Macrophages", "CDCs", "PDCs")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 9674321)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:45)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:45, seed.use = 632163)

# Save data
saveRDS(seuratObject, hc_colon_mnp_seurat_rds, compress = "gzip")

sessionInfo()