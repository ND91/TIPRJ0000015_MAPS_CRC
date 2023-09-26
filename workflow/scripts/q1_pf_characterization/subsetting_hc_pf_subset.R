#!/usr/bin/env Rscript
# The goal of this script is to filter the HC PF cells by the defined subset whereupon we perform clustering and dimension reduction using both gene expression and protein expression.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(Seurat)
library(SeuratDisk)
library(dplyr)

seurat_rds <- args[1]
subset <- args[2]
hc_pf_subset_seurat_rds <- args[3]

seuratObject <- readRDS(seurat_rds)

if(subset == "T"){
  seuratObject_HC_PF_subset <- seuratObject[,which(seuratObject@meta.data$manual_l1 %in% c("T"))]
  rna_pca_seed <- 67894123
  cite_pca_seed <- 51235
  rna_umap_seed <- 61263
  rna_npcs_test <- 100
  cite_dim <- 140
  rna_dim <- 32
} else if(subset == "NKILC"){
  seuratObject_HC_PF_subset <- seuratObject[,which(seuratObject@meta.data$manual_l1 %in% c("NK/ILC"))]
  rna_pca_seed <- 7894321
  cite_pca_seed <- 15215
  rna_umap_seed <- 143
  rna_npcs_test <- 100
  cite_dim <- 140
  rna_dim <- 43
} else if(subset == "B"){ # Few B cells only. Couldn't get this to work
  seuratObject_HC_PF_subset <- seuratObject[,which(seuratObject@meta.data$manual_l1 %in% c("B"))]
  rna_pca_seed <- 689743289
  cite_pca_seed <- 4321321
  rna_umap_seed <- 12351
  rna_npcs_test <- 76 
  cite_dim <- 38 
  rna_dim <- 45
} else if(subset == "Myeloid"){
  seuratObject_HC_PF_subset <- seuratObject[,which(seuratObject@meta.data$manual_l1 %in% c("Myeloid"))]
  rna_pca_seed <- 15233621
  cite_pca_seed <- 6321
  rna_umap_seed <- 61326312
  rna_npcs_test <- 100
  cite_dim <- 140
  rna_dim <- 54
} else if(subset == "MNP"){
  seuratObject_HC_PF_subset <- seuratObject[,which(seuratObject@meta.data$manual_l2 %in% c("Monocytes", "Macrophages", "CDCs", "PDCs"))]
  rna_pca_seed <- 789
  cite_pca_seed <- 67892143
  rna_umap_seed <- 632644
  rna_npcs_test <- 100
  cite_dim <- 140
  rna_dim <- 37
}

# Normalize, reduce, and recluster
seuratObject_HC_PF_subset_cleaned <- DietSeurat(seuratObject_HC_PF_subset, counts = T, data = T, scale.data = F)
seuratObject_HC_PF_subset_cleaned <- seuratObject_HC_PF_subset_cleaned[Matrix::rowSums(seuratObject_HC_PF_subset_cleaned) != 0, ]
seuratObject_HC_PF_subset_cleaned[["CITE"]] <- seuratObject_HC_PF_subset[["CITE"]]

## RNA
DefaultAssay(seuratObject_HC_PF_subset_cleaned) <- "RNA"
seuratObject_HC_PF_subset_cleaned <- SCTransform(seuratObject_HC_PF_subset_cleaned, conserve.memory = T)
seuratObject_HC_PF_subset_cleaned <- RunPCA(object = seuratObject_HC_PF_subset_cleaned, npcs = rna_npcs_test, seed.use = rna_pca_seed, reduction.name = "gex.pca")
seuratObject_HC_PF_subset_cleaned <- RunUMAP(seuratObject_HC_PF_subset_cleaned, reduction = 'gex.pca', dims = 1:rna_dim, assay = 'RNA', reduction.name = 'gex.umap', reduction.key = 'gexUMAP_')

## CITE
DefaultAssay(seuratObject_HC_PF_subset_cleaned) <- "CITE"
VariableFeatures(seuratObject_HC_PF_subset_cleaned) <- rownames(seuratObject_HC_PF_subset_cleaned[["CITE"]])
seuratObject_HC_PF_subset_cleaned <- NormalizeData(seuratObject_HC_PF_subset_cleaned, normalization.method = "CLR", assay = "CITE")
seuratObject_HC_PF_subset_cleaned <- ScaleData(seuratObject_HC_PF_subset_cleaned, assay = "CITE")
seuratObject_HC_PF_subset_cleaned <- RunPCA(object = seuratObject_HC_PF_subset_cleaned, npcs = cite_dim, seed.use = cite_pca_seed, reduction.name = "cite.pca", approx = FALSE)
seuratObject_HC_PF_subset_cleaned <- RunUMAP(seuratObject_HC_PF_subset_cleaned, reduction = 'cite.pca', dims = 1:cite_dim, assay = 'CITE', reduction.name = 'cite.umap', reduction.key = 'citeUMAP_')

seuratObject_HC_PF_subset_cleaned <- FindMultiModalNeighbors(seuratObject_HC_PF_subset_cleaned, reduction.list = list("gex.pca", "cite.pca"), dims.list = list(1:rna_dim, 1:cite_dim), modality.weight.name = "RNA.weight")
seuratObject_HC_PF_subset_cleaned <- FindClusters(seuratObject_HC_PF_subset_cleaned, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
seuratObject_HC_PF_subset_cleaned <- RunUMAP(seuratObject_HC_PF_subset_cleaned, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = rna_umap_seed)

seuratObject_HC_PF_subset_cleaned <- NormalizeData(seuratObject_HC_PF_subset_cleaned, assay = "CITE", normalization.method = "CLR")

# Save data
saveRDS(seuratObject_HC_PF_subset_cleaned, hc_pf_subset_seurat_rds, compress = "gzip")

sessionInfo()