#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 41.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(ggpubr)
require(ggplotify)

seurat_gcpmp_pf_rds <- args[1]
seurat_gcpmp_pf_mnp_rds <- args[2] 
seurat_gcpmp_pf_macrophages_rds <- args[3] #"resources/gastric_cancer/maps_gc_pf_macrophages_txpts.Rds"

# Immune

seurat_gcpmp_pf_tx <- readRDS(seurat_gcpmp_pf_rds)

seurat_gcpmp_pf

# MNP

# Macrophages

seurat_gcpmp_pf_macrophages <- readRDS(seurat_gcpmp_pf_macrophages_rds)

figR3P41 <- FeaturePlot(seurat_gcpmp_pf_macrophages, features = c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), label = T, ncol = 4)

ggsave(filename = "figR3P41.pdf", width = 18, height = 8, units = "in")
