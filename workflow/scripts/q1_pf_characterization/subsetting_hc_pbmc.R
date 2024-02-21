#!/usr/bin/env Rscript
# The goal of this script is to filter the HC PBMC whereupon we perform clustering and dimension reduction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
hc_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject_HC_PBMC <- seuratObject[,which(seuratObject@meta.data$Group %in% c("HC") & seuratObject@meta.data$Tissue %in% c("PBMC"))]

# Normalize, reduce, and recluster
seuratObject_HC_PBMC_cleaned <- SCTransform(seuratObject_HC_PBMC, verbose = FALSE, conserve.memory = TRUE)
seuratObject_HC_PBMC_cleaned <- RunPCA(object = seuratObject_HC_PBMC_cleaned, npcs = 100, seed.use = 843843)
seuratObject_HC_PBMC_cleaned <- FindNeighbors(seuratObject_HC_PBMC_cleaned, reduction = "pca", dims = 1:54)
seuratObject_HC_PBMC_cleaned <- FindClusters(seuratObject_HC_PBMC_cleaned, resolution = 0.5, verbose = FALSE)
seuratObject_HC_PBMC_cleaned <- RunUMAP(seuratObject_HC_PBMC_cleaned, dims = 1:54, seed.use = 71731727)

# Save data
saveRDS(seuratObject_HC_PBMC_cleaned, hc_seurat_rds_path, compress = "gzip")

sessionInfo()