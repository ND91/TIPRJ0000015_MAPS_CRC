#!/usr/bin/env Rscript

# The goal of this script is to combine the seuratObjects of PF (HC), liver, and colon.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))

seurat_pbmc_pf_rds <- args[1]
seurat_liver_rds <- args[2]
seurat_colon_rds <- args[3]
seurat_pf_liver_colon_rds <- args[4]

seurat_pbmc_pf <- readRDS(seurat_pbmc_pf_rds)
seurat_liver <- readRDS(seurat_liver_rds)
seurat_colon <- readRDS(seurat_colon_rds)

seurat_liver@meta.data <- seurat_liver@meta.data %>%
  dplyr::rename(manual_l0 = cell_type_MAPS_L0,
                manual_l1 = cell_type_MAPS_L1,
                manual_l2 = cell_type_MAPS_L2,
                manual_l3 = cell_type_MAPS_L3,
                Donor = patient) %>%
  dplyr::mutate(SampleID = paste0(Donor, "_Liver")) %>%
  dplyr::mutate(Tissue = "Liver",
                Section = "Liver",
                Study = "[liver study]",
                CellID = paste0(rownames(.), "_liver"))
seurat_liver <- RenameCells(seurat_liver, new.names = seurat_liver@meta.data$CellID)

seurat_colon@meta.data <- seurat_colon@meta.data %>%
  dplyr::rename(manual_l1 = cell_type_MAPS_L1,
                manual_l2 = cell_type_MAPS_L2,
                manual_l3 = cell_type_MAPS_L3,
                Donor = donor) %>%
  dplyr::mutate(Study = "[colon study]",
                Tissue = "Colon",
                Section = region,
                SampleID = paste0(Donor, "_", Tissue),
                CellID = paste0(rownames(.), "_colon"))
seurat_colon <- RenameCells(seurat_colon, new.names = seurat_colon@meta.data$CellID)

seurat_pbmc_pf@meta.data <- seurat_pbmc_pf@meta.data %>%
  dplyr::mutate(Study = "Current study",
                Section = Tissue)

seurat_pbmc_pf_liver_colon <- Reduce(merge, list(seurat_pbmc_pf, seurat_liver, seurat_colon))

# Normalize, reduce, and recluster
seurat_pbmc_pf_liver_colon <- DietSeurat(seurat_pbmc_pf_liver_colon, counts = T, data = T, scale.data = F)
seurat_pbmc_pf_liver_colon <- seurat_pbmc_pf_liver_colon[Matrix::rowSums(seurat_pbmc_pf_liver_colon) != 0, ]

## RNA
DefaultAssay(seurat_pbmc_pf_liver_colon) <- "RNA"
seurat_pbmc_pf_liver_colon <- SCTransform(seurat_pbmc_pf_liver_colon, conserve.memory = T)
seurat_pbmc_pf_liver_colon <- RunPCA(object = seurat_pbmc_pf_liver_colon, npcs = 100, seed.use = 679844132)
seurat_pbmc_pf_liver_colon <- FindNeighbors(seurat_pbmc_pf_liver_colon, reduction = "pca", dims = 1:73)
seurat_pbmc_pf_liver_colon <- FindClusters(seurat_pbmc_pf_liver_colon, resolution = 0.5, verbose = FALSE)
seurat_pbmc_pf_liver_colon <- RunUMAP(seurat_pbmc_pf_liver_colon, dims = 1:73, seed.use = 21346789)

# Save data
saveRDS(seurat_pbmc_pf_liver_colon, seurat_pf_liver_colon_rds, compress = "gzip")

sessionInfo()