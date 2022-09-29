#!/usr/bin/env Rscript
# The goal of this script is to select the PBMC-, PF-, and TX-derived cells from the donors that provided all three only, where we select for the CD45+, non-dead, non-proliferating, non-multiplet cells.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
pbmc_pf_tx_paired_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

discard_cellIDs <- unique(c(grep("multiplet", seuratObject@meta.data$manual_l1, ignore.case = T), 
                            grep("proliferating", seuratObject@meta.data$manual_l4, ignore.case = T), 
                            grep("dead/debris", seuratObject@meta.data$manual_l1, ignore.case = T),
                            grep("unknown", seuratObject@meta.data$manual_l1, ignore.case = T),
                            grep("endothelial", seuratObject@meta.data$manual_l1, ignore.case = T),
                            grep("epithelial", seuratObject@meta.data$manual_l1, ignore.case = T),
                            grep("mesenchymal", seuratObject@meta.data$manual_l1, ignore.case = T)
                         ))
seuratObject <- seuratObject[,-discard_cellIDs]

txdonors <- unique(seuratObject@meta.data$Donor[seuratObject@meta.data$Tissue == "TX"])

seuratObject <- seuratObject[,seuratObject@meta.data$Donor %in% txdonors]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 2349678)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:45)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunTSNE(seuratObject, dims = 1:45, seed.use = 8982364)

# Save data
saveRDS(seuratObject, pbmc_pf_tx_paired_seurat_rds_path, compress = "gzip")

sessionInfo()