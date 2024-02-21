#!/usr/bin/env Rscript

# This script will perform marker gene analyses in a one-vs-all comparison for all macrophages at manual_l4.

require(Seurat)
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

hc_pbmc_pf_seuratobject_rds <- args[1]
findmarkers_list_rds <- args[2]

seuratObject <- readRDS(hc_pbmc_pf_seuratobject_rds)
Idents(seuratObject) <- "manual_l4"

markers_list <- FindAllMarkers(seuratObject, only.pos = T)
markers_list <- split(markers_list, markers_list$cluster)

saveRDS(markers_list, findmarkers_list_rds, compress = "gzip")

sessionInfo()
