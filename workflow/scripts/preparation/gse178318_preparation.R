#!/usr/bin/env Rscript
# The goal of this script is to import, normalize, and (superficially) annotate the data from GSE178318

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

gse178318_matrix_mtx <- args[1]
gse178318_genes_tsv <- args[2]
gse178318_barcodes_tsv <- args[3]
gse178318_cellannotations_csv <- args[4]
seurat_rds <- args[5]

cnts <- ReadMtx(mtx = gse178318_matrix_mtx, 
                cells = gse178318_barcodes_tsv,
                features = gse178318_genes_tsv)

cellannotations <- read.csv(gse178318_cellannotations_csv)

seuratObject <- CreateSeuratObject(counts = cnts, min.cells = 3)

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(CellID = rownames(.),
                Donor = gsub("[ATCG]+_(.+)_.+$", "\\1", CellID),
                Tumor = gsub("[ATCG]+_.+_(.+)$", "\\1", CellID))

seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]
seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 678932417)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:77)
seuratObject <- FindClusters(seuratObject, resolution = 0.5, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:77, seed.use = 63122)
seuratObject[["percent_MT"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")

seuratObject@meta.data %>%
  dplyr::left_join(cellannotations, by = "CellID")

rownames(seuratObject@meta.data) <- seuratObject@meta.data$CellID

saveRDS(seuratObject, file = seurat_rds, compress = "gzip")

sessionInfo()