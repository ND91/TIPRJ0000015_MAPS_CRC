#!/usr/bin/env Rscript

# The goal of this script is to select all cells (not only immune) from PF derived from HC and CRC PM+ patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_crc_pmp_allcells_seurat_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Tissue %in% "PF",
                Group %in% c("CRC+", "HC"),
                !manual_l2 %in% c("Lineage negative HP+PRG4+", "Dead/debris", "Immune SMIM25+", "Lineage negative HP+PRG4+", "Multiplets", "NK/ILC proliferating", "T apoptotic", "T proliferating", "DC proliferating", "Endothelial", "HSPCs")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]
seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 52345)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:41)
seuratObject <- FindClusters(seuratObject, resolution = 1, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:41, seed.use = 32)

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l1l2 = manual_l1,
                manual_l1l2 = ifelse(manual_l2 %in% c("Myofibroblasts", "Fibroblasts"), "(Myo)Fibroblasts", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("Monocytes", "Macrophages"), "Mono-macs", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("CDCs", "PDCs"), "DCs", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("Granulocytes"), "Granulocytes", manual_l1l2))

# Save data
saveRDS(seuratObject, hc_crc_pmp_allcells_seurat_rds, compress = "gzip")

sessionInfo()