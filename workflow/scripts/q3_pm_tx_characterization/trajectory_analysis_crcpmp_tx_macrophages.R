#!/usr/bin/env Rscript

# This script will perform the TSCAN trajectory analysis on the monocytes and macrophages from PF and TX acquired from PM+ CRC (paired).

require(Seurat)
require(dplyr)
require(scater)
require(TSCAN)
require(tradeSeq)
require(slingshot)
require(scran)
require(scuttle)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
sce_rds <- args[2]
aggregated_lines_csv <- args[3]
tscan_rootmacrophagesvcan_rds <- args[4]

seuratObject <- readRDS(seurat_rds)
sce <- as.SingleCellExperiment(seuratObject)
colLabels(sce) <- colData(sce)$manual_l4

# Creating MST 
aggregated_clusters <- aggregateAcrossCells(sce, ids=colLabels(sce))
aggregated_centroids <- reducedDim(aggregated_clusters, "PCA")
aggregated_mst <- TSCAN::createClusterMST(aggregated_centroids, clusters = NULL)
aggregated_lines <- TSCAN::reportEdges(aggregated_clusters, mst=aggregated_mst, use.dimred="UMAP")
cellmapped_aggregatedmst <- mapCellsToEdges(sce, mst=aggregated_mst, use.dimred="PCA")

write.csv(aggregated_lines, aggregated_lines_csv)

# Rooting the MST
colData(sce)$entropy <- perCellEntropy(sce)

tscan_rootmacrophagesvcan <- orderCells(cellmapped_aggregatedmst, aggregated_mst, start="Macrophages VCAN+")
#colData(sce)$Pseudotime <- TSCAN::averagePseudotime(tscan_rootmacrophagesvcan) 

saveRDS(sce, sce_rds, compress = "gzip")
saveRDS(tscan_rootmacrophagesvcan, tscan_rootmacrophagesvcan_rds, compress = "gzip")