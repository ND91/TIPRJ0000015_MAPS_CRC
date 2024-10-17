#!/usr/bin/env Rscript

# The goal of this script is to select the TX-derived macrophages from the CRC patients.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
crc_tx_monocytes_macrophages_seurat_rds <- args[2]
crc_tx_monocytes_macrophages_cellmetadata_csv <- args[3]
crc_tx_monocytes_macrophages_counts_mtx <- args[4]
crc_tx_monocytes_macrophages_features_csv <- args[5]
crc_tx_monocytes_macrophages_umap_csv <- args[6]
crc_tx_monocytes_macrophages_pca_csv <- args[7]

seuratObject <- readRDS(seurat_rds)

txmonomac <- seuratObject@meta.data %>%
  dplyr::filter(Tissue == "TX",
                manual_l2 %in% c("Monocytes", "Macrophages")) %>%
  dplyr::pull(cellID)

seuratObject <- seuratObject[,seuratObject@meta.data$cellID %in% txmonomac]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]

seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 12363216)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:30)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:30, seed.use = 762589716)

# Save data
saveRDS(seuratObject, crc_tx_monocytes_macrophages_seurat_rds, compress = "gzip")

write.csv(seuratObject@meta.data, file = crc_tx_monocytes_macrophages_cellmetadata_csv, quote = F)
Matrix::writeMM(GetAssayData(seuratObject, assay = 'RNA', layer = 'counts'), file = crc_tx_monocytes_macrophages_counts_mtx)
write.table(data.frame('gene' = rownames(seuratObject)), file = crc_tx_monocytes_macrophages_features_csv, quote = F, row.names=F, col.names=F)
write.csv(Embeddings(seuratObject, reduction = "umap"), file = crc_tx_monocytes_macrophages_umap_csv, quote=F, row.names=F)
write.csv(Embeddings(seuratObject, reduction = "pca"), file = crc_tx_monocytes_macrophages_pca_csv, quote=F, row.names=F)

sessionInfo()