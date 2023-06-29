#!/usr/bin/env Rscript

# This script will perform the trajectory inference analysis on the monocytes and macrophages from PBMC, PF, and TX acquired from HC and PM+ CRC (paired).

require(Seurat)
require(dplyr)
require(scater)
require(tradeSeq)
require(slingshot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1]
sce_ss_rds <- args[2]

seuratObject <- readRDS(seurat_rds)

sce <- as.SingleCellExperiment(seuratObject, assay = "RNA")

sce <- slingshot(data = sce, 
                 reducedDim = 'UMAP', 
                 clusterLabels = "manual_l3",
                 start.clus = 'Classical monocytes', 
                 approx_points = 150)

saveRDS(sce, sce_ss_rds)
