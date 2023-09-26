#!/usr/bin/env Rscript
# The goal of this script is to append the manual curations.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(dplyr))

seurat_rds <- args[1]
curated_annotations_csv <- args[2]
seurat_curated_rds <- args[3]
seurat_metadata_csv <- args[4]

seuratObject <- readRDS(seurat_rds)
curated_annotations <- read.csv(curated_annotations_csv)

seuratObject_curated <- seuratObject

metadata <- seuratObject_curated@meta.data %>%
  dplyr::mutate(cellID = rownames(.)) %>%
  dplyr::left_join(curated_annotations, by = "cellID")

rownames(metadata) <- metadata$cellID

seuratObject_curated@meta.data <- metadata
seuratObject_curated[["CITE"]] <- seuratObject[["CITE"]]

# Save data
saveRDS(seuratObject_curated, seurat_curated_rds, compress = "gzip")
write.csv(seuratObject_curated@meta.data, seurat_metadata_csv)

sessionInfo()