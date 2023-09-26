#!/usr/bin/env Rscript
# The goal of this script is to convert the loom files into sce files and to rename the CellIDs to the correct format.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste0("Script needs at least 2 arguments. Current input is:", args))
}
loom <- args[1]
sce_rds <- args[2]

require(tidyverse)
require(LoomExperiment)
require(SingleCellExperiment)

scle <- import(loom, type="SingleCellLoomExperiment")
colData(scle) <- colData(scle) %>%
  data.frame() %>%
  mutate(CellID = paste0(gsub("-", "_", gsub("(.+):([ATCG]{16})x", "S\\1_\\2", CellID)), "-1")) %>%
  DataFrame()

sce <- as(scle, "SingleCellExperiment")

saveRDS(sce, sce_rds, compress = "gzip")