#!/usr/bin/env Rscript
# The goal of this script is to merge all the datasets into one large Seurat object.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste0("Script needs at least 5 arguments. Current input is:", args))
}
seurat_merged_rds <- args[1]

seurat_paths <- args[2:length(args)] 

t0 <- Sys.time()

seurat_list <- lapply(seurat_paths, readRDS)

seurat_merged <- Reduce(merge, seurat_list)

paste0("Reading and merging all seurat objects took ", as.integer(as.numeric(Sys.time() - t0, unit = "secs")), " seconds")

DefaultAssay(seurat_merged) <- "RNA"

saveRDS(seurat_merged, seurat_merged_rds, compress = "gzip")

sessionInfo()

