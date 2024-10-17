#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 2, point 1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(pheatmap)
require(ggpubr)

seurat_rds <- args[1] #"output/subsets/live_singlet_nonproliferating_SeuratObject.Rds"
macrophages_pf_tx_seurat_rds <- args[2] #"output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_SeuratObject.Rds"

seuratObject <- readRDS(seurat_rds)

selected_cells <- seuratObject@meta.data %>%
  dplyr::filter(Tissue %in% c("PF", "TX"),
  )

pf_tx_seurat <- seuratObject[,which(seuratObject@meta.data$CellID %in% selected_cells$CellID)]

Idents(pf_tx_seurat) <- "manual_l2"

figR2P1A <- VlnPlot(pf_tx_seurat, "VSIG4", group.by = "manual_l2")
figR2P1B <- FeaturePlot(pf_tx_seurat, "VSIG4", split.by = "Tissue", label = T)

macrophage_cells <- pf_tx_seurat@meta.data %>%
  dplyr::filter(manual_l2 %in% c("Macrophages"),
  )

macrophages_pf_tx_seurat <- readRDS(macrophages_pf_tx_seurat_rds)

figR2P1C <- DimPlot(macrophages_pf_tx_seurat, split.by = "Tissue", label = T, group.by = "manual_l3")
figR2P1D <- VlnPlot(macrophages_pf_tx_seurat, "VSIG4", group.by = "manual_l3") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

Idents(macrophages_pf_tx_seurat) <- "manual_l3"

figR2P1E <- FeaturePlot(macrophages_pf_tx_seurat, "VSIG4", split.by = "Tissue", label = T)

figR2P1AB <- ggarrange(figR2P1A, figR2P1B, nrow = 1, ncol = 2, labels = c("A", "B"))
figR2P1CDE <- ggarrange(figR2P1C, figR2P1D, figR2P1E[[1]], figR2P1E[[2]], nrow = 1, ncol = 4, labels = c("C", "D", "E", ""), widths = c(1.5, 0.5, 0.5, 0.5))

figR2P1 <- ggarrange(figR2P1AB, figR2P1CDE, nrow = 2, ncol = 1)

ggsave(filename = "figR2P1.pdf", width = 22.5, height = 10, units = "in")
