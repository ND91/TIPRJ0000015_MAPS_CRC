#!/usr/bin/env Rscript

# The goal of this script is to select the TX-derived immune cells from the CRC patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crc_tx_immune_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

txcells <- seuratObject@meta.data %>%
  dplyr::filter(manual_l1 %in% c("B", "Myeloid", "NK/ILC", "T"),
                #!manual_l2 %in% c("Granulocytes"),
                !manual_l4 %in% c("T apoptotic", "Platelets", "NK/ILC proliferating", "T proliferating", "CD8 T memory proliferating"),
                Tissue %in% c("TX")) %>%
  dplyr::pull(cellID)

seuratObject <- seuratObject[,seuratObject@meta.data$cellID %in% txcells]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 1432)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:68)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:68, seed.use = 631623)

# Save data
saveRDS(seuratObject, crc_tx_immune_seurat_rds, compress = "gzip")

sessionInfo()