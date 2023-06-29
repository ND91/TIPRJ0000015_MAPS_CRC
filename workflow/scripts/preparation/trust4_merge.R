#!/usr/bin/env Rscript
# The goal of this script is to merge all the AIRR files into one large AIRR CSV and annotate it with the cell metadata as obtained from the seuratobject.

library(dplyr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste0("Script needs at least 5 arguments. Current input is:", args))
}
airr_merged_csv <- args[1]
seuratObject_rds <- args[2]

seuratObject <- readRDS(seuratObject_rds)

airr_annotated_csvs <- args[3:length(args)] 

airr_list <- lapply(airr_annotated_csvs, read.csv)

airr_merged <- do.call(rbind, airr_list)

airr_merged <- airr_merged %>%
  dplyr::inner_join(seuratObject@meta.data, by = "cellID")

write.csv(airr_merged, airr_merged_csv, row.names = F)

sessionInfo()