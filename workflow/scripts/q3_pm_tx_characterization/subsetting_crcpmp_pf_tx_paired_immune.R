#!/usr/bin/env Rscript

# The goal of this script is to select the PF- and TX-derived cells from the CRC patients that provided all three samples.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crc_pf_tx_paired_immune_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

txdonors <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX") %>%
  dplyr::pull(Donor) %>%
  unique()

seuratObject <- seuratObject[,seuratObject@meta.data$Donor %in% txdonors]

selected_cells <- seuratObject@meta.data %>%
  dplyr::filter(manual_l1 %in% c("B", "Myeloid", "NK/ILC", "T"),
                #!manual_l2 %in% c("Granulocytes"),
                !manual_l4 %in% c("T apoptotic", "Platelets", "NK/ILC proliferating", "T proliferating", "CD8 T memory proliferating"),
                Tissue %in% c("PF", "TX")
  )

seuratObject_filtered <- seuratObject[,seuratObject@meta.data$cellID %in% selected_cells$cellID]

# Normalize, reduce, and recluster
seuratObject_filtered <- DietSeurat(seuratObject_filtered, counts = T, data = T, scale.data = F)
seuratObject_filtered <- seuratObject_filtered[Matrix::rowSums(seuratObject_filtered) != 0, ]

seuratObject_filtered <- SCTransform(seuratObject_filtered, conserve.memory = T)
seuratObject_filtered <- RunPCA(object = seuratObject_filtered, npcs = 100, seed.use = 631271)
seuratObject_filtered <- FindNeighbors(seuratObject_filtered, reduction = "pca", dims = 1:62)
seuratObject_filtered <- FindClusters(seuratObject_filtered, resolution = 0.5, verbose = FALSE)
seuratObject_filtered <- RunUMAP(seuratObject_filtered, dims = 1:62, seed.use = 7127333)

# Save data
saveRDS(seuratObject_filtered, crc_pf_tx_paired_immune_seurat_rds, compress = "gzip")

sessionInfo()