#!/usr/bin/env Rscript

# The goal of this script is to select the PBMC and PF monocytes and macrophages from paired CRC PM+ patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

txdonor <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX") %>%
  dplyr::pull(Donor) %>%
  unique()

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Donor %in% c(txdonor),
                manual_l2 %in% c("Monocytes", "Macrophages"),
                !manual_l4 %in% c("NK/ILC proliferating", "T proliferating", "DC proliferating", "Macrophages proliferating", "CD8 T memory proliferating", "T apoptotic", "Platelets")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 39780)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:48)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:48, seed.use = 361236)

# Save data
saveRDS(seuratObject, crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_rds, compress = "gzip")

sessionInfo()