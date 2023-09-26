#!/usr/bin/env Rscript
# The goal of this script is to import and reannotate the loom files using the annotations provided for the SeuratObject.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

loom <- args[1]
seuratObject_rds <- args[2]
loom_annotated <- args[3]

require(tidyverse)
require(LoomExperiment)
require(SingleCellExperiment)
require(Seurat)

scle <- LoomExperiment::import(loom, type="SingleCellLoomExperiment")

seuratObject <- readRDS(seuratObject_rds)

colData(scle) <- colData(scle) %>%
  data.frame() %>%
  mutate(CellID_annotated = paste0(gsub("-", "_", gsub("(.+):([ATCG]{16})x", "S\\1_\\2", CellID)), "-1")) %>%
  DataFrame()

# colnames(scle) <- colData(scle)$CellID
# rownames(scle) <- rowData(scle)$Accession

scle_selected <- scle[,which(colData(scle)$CellID %in% seuratObject$CellID)]

colData(scle_selected) <- colData(scle_selected) %>%
  data.frame() %>%
  dplyr::left_join(seuratObject@meta.data, by = c("CellID_annotated" = "CellID")) %>%
  DataFrame()

LoomExperiment::export(object = scle_selected, con = loom_annotated)

sessionInfo()