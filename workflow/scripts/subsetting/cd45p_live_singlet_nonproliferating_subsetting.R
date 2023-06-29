#!/usr/bin/env Rscript
# The goal of this script is to filter the cleaned (CD45+, live, singlet and non-proliferating) cells whereupon we perform clustering and dimension reduction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
cleaned_seurat_rds_path <- args[2]

seuratObject <- readRDS(seurat_rds_path)

selected_cells <- seuratObject@meta.data %>%
  dplyr::filter(manual_l1 %in% c("B", "Myeloid", "NK/ILC", "T"),
                !manual_l2 %in% c("Granulocytes"),
                !manual_l4 %in% c("NK/ILC proliferating", "T proliferating", "DC proliferating", "Macrophages proliferating", "CD8 T memory proliferating", "T apoptotic", "Platelets"),
  )

seuratObject_filtered <- seuratObject[,seuratObject@meta.data$cellID %in% selected_cells$cellID]

# Normalize, reduce, and recluster
seuratObject_filtered_cleaned <- DietSeurat(seuratObject_filtered, counts = T, data = T, scale.data = F)
seuratObject_filtered_cleaned <- seuratObject_filtered_cleaned[Matrix::rowSums(seuratObject_filtered_cleaned) != 0, ]

seuratObject_filtered_cleaned <- SCTransform(seuratObject_filtered_cleaned, conserve.memory = T)
seuratObject_filtered_cleaned <- RunPCA(object = seuratObject_filtered_cleaned, npcs = 100, seed.use = 41326789)
seuratObject_filtered_cleaned <- FindNeighbors(seuratObject_filtered_cleaned, reduction = "pca", dims = 1:67)
seuratObject_filtered_cleaned <- FindClusters(seuratObject_filtered_cleaned, resolution = 0.5, verbose = FALSE)
seuratObject_filtered_cleaned <- RunUMAP(seuratObject_filtered_cleaned, dims = 1:67, seed.use = 6712349)
seuratObject_filtered_cleaned[["CITE"]] <- seuratObject_filtered[["CITE"]]

# Save data
saveRDS(seuratObject_filtered_cleaned, cleaned_seurat_rds_path, compress = "gzip")

sessionInfo()