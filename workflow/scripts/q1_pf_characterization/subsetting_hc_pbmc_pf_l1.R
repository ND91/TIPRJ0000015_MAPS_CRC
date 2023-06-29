#!/usr/bin/env Rscript
# The goal of this script is to filter the HC PBMC and PF for the l1 type whereupon we perform clustering and dimension reduction using both gene expression and protein expression.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(Seurat)
library(SeuratDisk)
library(dplyr)

seurat_rds <- args[1]
l1 <- args[2]
hc_pbmc_pf_l1_seurat_rds <- args[3]

seuratObject <- readRDS(seurat_rds)

cluster_params <- data.frame(
  ndims = dplyr::case_when(
    l1 == "Myeloid" ~ 43,
    l1 == "B" ~ 32),
  pcaseed = dplyr::case_when(
    l1 == "Myeloid" ~ 26978,
    l1 == "B" ~ 1536928),
  umapseed = dplyr::case_when(
    l1 == "Myeloid" ~ 42,
    l1 == "B" ~ 416923))

seuratObject_subset <- seuratObject[,which(seuratObject@meta.data$manual_l1 == l1)]

# Normalize, reduce, and recluster
seuratObject_subset_cleaned <- DietSeurat(seuratObject_subset, counts = T, data = T, scale.data = F)
seuratObject_subset_cleaned <- seuratObject_subset_cleaned[Matrix::rowSums(seuratObject_subset_cleaned) != 0, ]
seuratObject_subset_cleaned[["CITE"]] <- seuratObject_subset[["CITE"]]

## RNA
DefaultAssay(seuratObject_subset_cleaned) <- "RNA"
seuratObject_subset_cleaned <- SCTransform(seuratObject_subset_cleaned, conserve.memory = T)
seuratObject_subset_cleaned <- RunPCA(object = seuratObject_subset_cleaned, npcs = 100, seed.use = cluster_params$pcaseed[1])
seuratObject_subset_cleaned <- FindNeighbors(seuratObject_subset_cleaned, reduction = "pca", dims = 1:cluster_params$ndims[1])
seuratObject_subset_cleaned <- FindClusters(seuratObject_subset_cleaned, resolution = 0.5, verbose = FALSE)
seuratObject_subset_cleaned <- RunUMAP(seuratObject_subset_cleaned, dims = 1:cluster_params$ndims[1], seed.use = cluster_params$umapseed[1])

# Save data
saveRDS(seuratObject_subset_cleaned, hc_pbmc_pf_B_seurat_rds_path, compress = "gzip")

sessionInfo()