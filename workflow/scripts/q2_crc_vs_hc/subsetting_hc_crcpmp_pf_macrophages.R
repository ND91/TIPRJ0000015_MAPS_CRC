#!/usr/bin/env Rscript

# The goal of this script is to select the PF macrophages from HC and CRC PM+ patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_crcpmp_pf_macrophages_seurat_rds <- args[2]
hc_crcpmp_pf_macrophages_cellmetadata_csv <- args[3]
hc_crcpmp_pf_macrophages_counts_mtx <- args[4]
hc_crcpmp_pf_macrophages_features_csv <- args[5]
hc_crcpmp_pf_macrophages_umap_csv <- args[6]
hc_crcpmp_pf_macrophages_pca_csv <- args[7]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Group %in% c("CRC+", "HC"),
                Tissue %in% c("PF"),
                manual_l2 %in% c("Macrophages")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 4132)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:47)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:47, seed.use = 5132)

# Save data
saveRDS(seuratObject, hc_crcpmp_pf_macrophages_seurat_rds, compress = "gzip")

write.csv(seuratObject@meta.data, file = hc_crcpmp_pf_macrophages_cellmetadata_csv, quote = F)
Matrix::writeMM(GetAssayData(seuratObject, assay = 'RNA', layer = 'counts'), file = hc_crcpmp_pf_macrophages_counts_mtx)
write.table(data.frame('gene' = rownames(seuratObject)), file = hc_crcpmp_pf_macrophages_features_csv, quote = F, row.names=F, col.names=F)
write.csv(Embeddings(seuratObject, reduction = "umap"), file = hc_crcpmp_pf_macrophages_umap_csv, quote=F, row.names=F)
write.csv(Embeddings(seuratObject, reduction = "pca"), file = hc_crcpmp_pf_macrophages_pca_csv, quote=F, row.names=F)

sessionInfo()