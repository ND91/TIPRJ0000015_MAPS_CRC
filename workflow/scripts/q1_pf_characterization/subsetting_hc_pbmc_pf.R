#!/usr/bin/env Rscript
# The goal of this script is to filter the HC PBMC and PF whereupon we perform clustering and dimension reduction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Harmony))

seurat_rds_path <- args[1]
hc_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject_HC <- seuratObject[,which(seuratObject@meta.data$Group %in% c("HC"))]

# Normalize, reduce, and recluster
seuratObject_HC_cleaned <- SCTransform(seuratObject_HC, verbose = FALSE, conserve.memory = TRUE)
seuratObject_HC_cleaned <- RunPCA(object = seuratObject_HC_cleaned, npcs = 100, seed.use = 2349178)
seuratObject_HC_cleaned <- FindNeighbors(seuratObject_HC_cleaned, reduction = "pca", dims = 1:43)
seuratObject_HC_cleaned <- FindClusters(seuratObject_HC_cleaned, resolution = 0.5, verbose = FALSE)
# seuratObject_HC_cleaned <- RunUMAP(seuratObject_HC_cleaned, dims = 1:67, seed.use = 21346789)
seuratObject_HC_cleaned[["CITE"]] <- seuratObject_HC[["CITE"]]

seuratObject_HC_cleaned <- NormalizeData(seuratObject_HC_cleaned, assay = "CITE", normalization.method = "CLR")

DefaultAssay(seuratObject_HC_cleaned) <- "RNA"

# Harmony batch effect correction for tissue
seuratObject_HC_cleaned <- RunHarmony(seuratObject_HC_cleaned, "Tissue")
seuratObject_HC_cleaned <- RunUMAP(seuratObject_HC_cleaned, dims = 1:45, seed.use = 398772, reduction = "harmony")

# Save data
saveRDS(seuratObject_HC_cleaned, hc_seurat_rds_path, compress = "gzip")

sessionInfo()