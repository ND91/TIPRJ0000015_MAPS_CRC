#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrastr))
suppressPackageStartupMessages(require(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
umap_macrophages_manuall3 <- args[2]
umap_macrophages_markers <- args[3]
boxplot_macrophages_markers <- args[4]

seuratObject <- readRDS(seurat_rds_path)

macrophages_seuratObject <- seuratObject[,seuratObject@meta.data$manual_l2 == "Macrophages"]
macrophages_seuratObject_cleaned <- SCTransform(macrophages_seuratObject, verbose = FALSE, conserve.memory = TRUE)
macrophages_seuratObject_cleaned <- RunPCA(object = macrophages_seuratObject_cleaned, npcs = 100, seed.use = 1437280)
macrophages_seuratObject_cleaned <- FindNeighbors(macrophages_seuratObject_cleaned, reduction = "pca", dims = 1:47)
macrophages_seuratObject_cleaned <- FindClusters(macrophages_seuratObject_cleaned, resolution = 0.5, verbose = FALSE)
macrophages_seuratObject_cleaned <- RunUMAP(macrophages_seuratObject_cleaned, dims = 1:47, seed.use = 697815)
macrophages_seuratObject_cleaned[["CITE"]] <- macrophages_seuratObject[["CITE"]]

sessionInfo()