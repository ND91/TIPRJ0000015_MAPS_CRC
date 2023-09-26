#!/usr/bin/env Rscript
# The goal of this script is to filter the HC PBMC and PF whereupon we perform clustering and dimension reduction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
hc_pf_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

seuratObject_HC_PF <- seuratObject[,which(seuratObject@meta.data$Group %in% c("HC") & seuratObject@meta.data$Tissue %in% c("PF"))]

# Normalize, reduce, and recluster
seuratObject_HC_PF_cleaned <- DietSeurat(seuratObject_HC_PF, counts = T, data = T, scale.data = F)
seuratObject_HC_PF_cleaned <- seuratObject_HC_PF_cleaned[Matrix::rowSums(seuratObject_HC_PF_cleaned) != 0, ]
seuratObject_HC_PF_cleaned[["CITE"]] <- seuratObject_HC_PF[["CITE"]]

## RNA
DefaultAssay(seuratObject_HC_PF_cleaned) <- "RNA"
seuratObject_HC_PF_cleaned <- SCTransform(seuratObject_HC_PF_cleaned, conserve.memory = T)
seuratObject_HC_PF_cleaned <- RunPCA(object = seuratObject_HC_PF_cleaned, npcs = 100, seed.use = 2349178, reduction.name = "gex.pca")

## CITE
DefaultAssay(seuratObject_HC_PF_cleaned) <- "CITE"
seuratObject_HC_PF_cleaned <- SCTransform(seuratObject_HC_PF_cleaned, conserve.memory = T)
seuratObject_HC_PF_cleaned <- RunPCA(object = seuratObject_HC_PF_cleaned, npcs = 140, seed.use = 16794832, reduction.name = "cite.pca")

seuratObject_HC_PF_cleaned <- FindMultiModalNeighbors(seuratObject_HC_PF_cleaned, reduction.list = list("gex.pca", "cite.pca"), dims.list = list(1:67, 1:140), modality.weight.name = "RNA.weight")
seuratObject_HC_PF_cleaned <- FindClusters(seuratObject_HC_PF_cleaned, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
seuratObject_HC_PF_cleaned <- RunUMAP(seuratObject_HC_PF_cleaned, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 67148392)
seuratObject_HC_PF_cleaned <- RunUMAP(seuratObject_HC_PF_cleaned, reduction = 'gex.pca', dims = 1:67, assay = 'RNA', reduction.name = 'gex.umap', reduction.key = 'gexUMAP_')
seuratObject_HC_PF_cleaned <- RunUMAP(seuratObject_HC_PF_cleaned, reduction = 'cite.pca', dims = 1:140, assay = 'CITE', reduction.name = 'cite.umap', reduction.key = 'citeUMAP_')

seuratObject_HC_PF_cleaned <- NormalizeData(seuratObject_HC_PF_cleaned, assay = "CITE", normalization.method = "CLR")

# Save data
saveRDS(seuratObject_HC_PF_cleaned, hc_pf_seurat_rds_path, compress = "gzip")

sessionInfo()