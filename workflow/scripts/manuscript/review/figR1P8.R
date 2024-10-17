#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 8.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(ggpubr)

seurat_rds <- args[1] #"output/q3_pm_tx_characterization/subsets/hc_crcpmp_pf_tx_macrophages_SeuratObject.Rds"

seuratObject <- readRDS(seurat_rds)

seuratObject@meta.data$Tissue_group <- factor(paste0(seuratObject@meta.data$Tissue, " ", seuratObject@meta.data$Group), levels = c("PF HC", "PF CRC+", "TX CRC+"))

figR1P8A <- DimPlot(seuratObject, group.by = "Tissue_group")
figR1P8B <- DimPlot(seuratObject, group.by = "Tissue", split.by = "Group")

figR1P8 <- ggarrange(figR1P8A, figR1P8B, nrow = 1, ncol = 2, widths = c(1.25, 2))

ggsave(filename = "figR1P8.pdf", width = 18, height = 6, units = "in")
