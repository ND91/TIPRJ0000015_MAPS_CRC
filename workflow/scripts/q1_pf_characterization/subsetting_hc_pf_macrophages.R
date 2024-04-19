#!/usr/bin/env Rscript

# The goal of this script is to select the PF macrophages from HC patients and store it both as Seurat object in RDS format and separate data tables in CSV format.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)

seurat_rds <- args[1]
hc_pf_macrophages_seurat_rds <- args[2]
hc_pf_macrophages_cellmetadata_csv <- args[3]
hc_pf_macrophages_counts_mtx <- args[4]
hc_pf_macrophages_features_csv <- args[5]
hc_pf_macrophages_umap_csv <- args[6]
hc_pf_macrophages_pca_csv <- args[7]

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(manual_l2 %in% c("Macrophages")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]
seuratObject_cleaned <- seuratObject

# Normalize, reduce, and recluster
seuratObject_cleaned <- DietSeurat(seuratObject_cleaned, counts = T, data = T, scale.data = F)
seuratObject_cleaned <- seuratObject_cleaned[Matrix::rowSums(seuratObject_cleaned) != 0, ]
seuratObject_cleaned[["CITE"]] <- seuratObject[["CITE"]]

## RNA
DefaultAssay(seuratObject_cleaned) <- "RNA"
seuratObject_cleaned <- SCTransform(seuratObject_cleaned, conserve.memory = T)
seuratObject_cleaned <- RunPCA(object = seuratObject_cleaned, npcs = 100, seed.use = 67984321, reduction.name = "gex.pca")

## CITE
DefaultAssay(seuratObject_cleaned) <- "CITE"
seuratObject_cleaned <- SCTransform(seuratObject_cleaned, conserve.memory = T)
seuratObject_cleaned <- RunPCA(object = seuratObject_cleaned, npcs = 140, seed.use = 7869432, reduction.name = "cite.pca")

seuratObject_cleaned <- FindMultiModalNeighbors(seuratObject_cleaned, reduction.list = list("gex.pca", "cite.pca"), dims.list = list(1:59, 1:140), modality.weight.name = "RNA.weight")
seuratObject_cleaned <- FindClusters(seuratObject_cleaned, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
seuratObject_cleaned <- RunUMAP(seuratObject_cleaned, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 67148392)
seuratObject_cleaned <- RunUMAP(seuratObject_cleaned, reduction = 'gex.pca', dims = 1:59, assay = 'RNA', reduction.name = 'gex.umap', reduction.key = 'gexUMAP_')
seuratObject_cleaned <- RunUMAP(seuratObject_cleaned, reduction = 'cite.pca', dims = 1:140, assay = 'CITE', reduction.name = 'cite.umap', reduction.key = 'citeUMAP_')

seuratObject_cleaned <- NormalizeData(seuratObject_cleaned, assay = "CITE", normalization.method = "CLR")

# Save data
saveRDS(seuratObject_cleaned, hc_pf_macrophages_seurat_rds, compress = "gzip")

write.csv(seuratObject_cleaned@meta.data, file = hc_pf_macrophages_cellmetadata_csv, quote = F)
Matrix::writeMM(GetAssayData(seuratObject_cleaned, assay = 'RNA', layer = 'counts'), file = hc_pf_macrophages_counts_mtx)
write.table(data.frame('gene' = rownames(seuratObject_cleaned)), file = hc_pf_macrophages_features_csv, quote = F, row.names=F, col.names=F)
write.csv(Embeddings(seuratObject_cleaned, reduction = "gex.umap"), file = hc_pf_macrophages_umap_csv, quote=F, row.names=F)
write.csv(Embeddings(seuratObject_cleaned, reduction = "gex.pca"), file = hc_pf_macrophages_pca_csv, quote=F, row.names=F)

sessionInfo()