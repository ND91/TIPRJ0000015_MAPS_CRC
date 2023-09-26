#!/usr/bin/env Rscript

# The goal of this script is to combine the seuratObjects of PF macrophages (PM-CRC), with macrophages from colon and liver (CRC).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_crc_pf_macrophages_rds <- args[1]
seurat_crc_liver_colon_rds <- args[2]
seurat_crc_pf_liver_colon_rds <- args[3]

seurat_crc_pf_macrophages <- readRDS(seurat_crc_pf_macrophages_rds)
seurat_crc_liver_colon <- readRDS(seurat_crc_liver_colon_rds)

seurat_crc_liver_colon@meta.data <- seurat_crc_liver_colon@meta.data %>%
  dplyr::mutate(Tissue = case_when(
    Tumor == "CRC" ~ "Colon", 
    Tumor == "LM" ~ "Liver",
    Tumor == "PBMC" ~ "PBMC"),
    Study = "GSE178318") %>%
  dplyr::select(-Tumor)

seurat_crc_liver_colon_macrophages <- seurat_crc_liver_colon[,seurat_crc_liver_colon@meta.data$manual_l2 == "Macrophages" & seurat_crc_liver_colon@meta.data$Tissue %in% c("Colon", "Liver")]

seurat_crc_pf_macrophages@meta.data <- seurat_crc_pf_macrophages@meta.data %>%
  dplyr::mutate(Study = "Current study")

seurat_pf_liver_colon_macrophages <- merge(seurat_crc_pf_macrophages, seurat_crc_liver_colon_macrophages)

# Normalize, reduce, and recluster
seurat_pf_liver_colon_macrophages <- DietSeurat(seurat_pf_liver_colon_macrophages, counts = T, data = T, scale.data = F)
seurat_pf_liver_colon_macrophages <- seurat_pf_liver_colon_macrophages[Matrix::rowSums(seurat_pf_liver_colon_macrophages) != 0, ]

## RNA
DefaultAssay(seurat_pf_liver_colon_macrophages) <- "RNA"
seurat_pf_liver_colon_macrophages <- SCTransform(seurat_pf_liver_colon_macrophages, conserve.memory = T)
seurat_pf_liver_colon_macrophages <- RunPCA(object = seurat_pf_liver_colon_macrophages, npcs = 100, seed.use = 5978243)
seurat_pf_liver_colon_macrophages <- FindNeighbors(seurat_pf_liver_colon_macrophages, reduction = "pca", dims = 1:42)
seurat_pf_liver_colon_macrophages <- FindClusters(seurat_pf_liver_colon_macrophages, resolution = 0.5, verbose = FALSE)
seurat_pf_liver_colon_macrophages <- RunUMAP(seurat_pf_liver_colon_macrophages, dims = 1:42, seed.use = 62346)

# Save data
saveRDS(seurat_pf_liver_colon_macrophages, seurat_crc_pf_liver_colon_rds, compress = "gzip")

sessionInfo()