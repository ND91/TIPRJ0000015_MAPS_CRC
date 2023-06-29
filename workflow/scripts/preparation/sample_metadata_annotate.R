#!/usr/bin/env Rscript
# The goal of this script is to append the patient metadata.

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
sample_metadata_xlsx <- args[2]
seurat_sample_metadata_annotated_rds <- args[3]

seuratObject <- readRDS(seurat_rds)

sample_metadata <- readxl::read_excel(sample_metadata_xlsx)

seuratObject_annotated <- seuratObject

metadata <- seuratObject@meta.data %>%
  dplyr::left_join(sample_metadata, by = "RunID")

metadata <- metadata[,order(colnames(metadata))]

rownames(metadata) <- colnames(seuratObject_annotated)

seuratObject_annotated@meta.data <- metadata

if("ADT" %in% Assays(seuratObject)){
  seuratObject_annotated[["CITE"]] <- seuratObject[["CITE"]]
}

DefaultAssay(seuratObject_annotated) <- "RNA"

saveRDS(seuratObject_annotated, seurat_sample_metadata_annotated_rds, compress = "gzip")

sessionInfo()
