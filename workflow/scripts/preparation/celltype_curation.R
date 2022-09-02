#!/usr/bin/env Rscript
# The goal of this script is to append the manual curations.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds_path <- args[1]
curated_annotations_csv_path <- args[2]
seurat_curated_rds_path <- args[3]

seuratObject <- readRDS(seurat_rds_path)
curated_annotations <- read.csv(curated_annotations_csv_path)[,-1]

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(cellID = rownames(.)) %>%
  dplyr::left_join(curated_annotations, by = "cellID")

rownames(seuratObject@meta.data) <- seuratObject@meta.data$cellID

# Save data
saveRDS(seuratObject, seurat_curated_rds_path, compress = "gzip")

sessionInfo()