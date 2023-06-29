#!/usr/bin/env Rscript
# The goal of this script is to perform clustering and dimension reduction on the non-proliferating myeloid samples from the myeloid fraction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
pf_atlas_path <- args[2]
pbmc_pf_tx_seurat_rds_path <- args[3]

seuratObject <- readRDS(seurat_rds_path)
seuratObject <- seuratObject[,seuratObject@meta.data$manual_l1 == "Myeloid" & seuratObject@meta.data$Tissue %in% c("PF", "PBMC")]

pf_atlas_seuratObject <- readRDS(pf_atlas_path)
pf_atlas_seuratObject@meta.data <- pf_atlas_seuratObject@meta.data %>%
  dplyr::rename(manual_l1 = celltype_l1,
                manual_l2 = celltype_l2,
                manual_l3 = celltype_l3,
                manual_l4 = celltype_l4) %>%
  dplyr::mutate(Group = "HC")

pf_atlas_seuratObject <- pf_atlas_seuratObject[,pf_atlas_seuratObject@meta.data$manual_l1 == "Myeloid"]

# Merge

seuratObject <- merge(seuratObject, pf_atlas_seuratObject)

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 967831)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:34)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunTSNE(seuratObject, dims = 1:34, seed.use = 532151)

# Save data
saveRDS(seuratObject, pbmc_pf_tx_seurat_rds_path, compress = "gzip")

sessionInfo()