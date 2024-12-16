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
hc_pf_seurat_rds <- args[2] #"output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds"
macrophages_pf_tx_seurat_rds <- args[3] #"output/q3_pm_tx_characterization/subsets/crcpmp_pf_tx_paired_macrophages_SeuratObject.Rds"
celltype_markers_xlsx <- args[4] #"config/order/celltype_markers.xlsx"

hc_pf_seurat <- readRDS(hc_pf_seurat_rds)

celltype_order_pf_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(hc_pf_seurat@meta.data[,"manual_l3"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l3_colors_pf <- celltype_order_pf_l3$color
names(manual_l3_colors_pf) <- celltype_order_pf_l3$celltype

# PF and TX

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

# HC

violinplot_hc_pf_vsig4_ggplotobj <- 
  data.frame(CB = colnames(hc_pf_seurat),
             Embeddings(hc_pf_seurat[["wnn.umap"]]),
             hc_pf_seurat@meta.data,
             expr = GetAssayData(hc_pf_seurat, assay = "RNA")["VSIG4",],
             Feature = "VSIG4") %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset)) %>%
  dplyr::filter(manual_l2 %in% "Macrophages") %>%
  ggplot(aes(x = celltype, y = expr, fill = celltype)) +
  geom_violin(trim = T, scale = "width", bounds = c(0, Inf)) +
  geom_jitter(size = 0.01) +
  labs(y = "nUMIs") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_fill_manual(values = manual_l3_colors_pf) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "lines"))

# ggsave(filename = "violinplot_hc_pf_vsig4.pdf", width = 10, height = 4, units = "in")
ggsave(filename = "violinplot_hc_pf_macrophages_vsig4.pdf", width = 3, height = 4, units = "in")

# PM-CRC PF and PM macrophages

macrophages_pf_tx_seurat <- readRDS(macrophages_pf_tx_seurat_rds)

figR2P1C <- DimPlot(macrophages_pf_tx_seurat, split.by = "Tissue", label = T, group.by = "manual_l3")
figR2P1D <- VlnPlot(macrophages_pf_tx_seurat, "VSIG4", group.by = "manual_l3") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

Idents(macrophages_pf_tx_seurat) <- "manual_l3"


# Percentage VSIG4-expressing macrophages

data.frame(hc_pf_seurat@meta.data, 
           VSIG4_expressing = GetAssayData(hc_pf_seurat, assay = "RNA")["VSIG4",]>0) %>%
  dplyr::filter(manual_l2 == "Macrophages") %>%
  dplyr::group_by(Donor, manual_l3, VSIG4_expressing) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Donor, manual_l3) %>%
  dplyr::mutate(percentage = n/(sum(n))*100) %>%
  dplyr::filter(VSIG4_expressing == T) %>%
  readr::write_csv("percentage_vsig4expressing_hc_pf_macrophages.csv")

figR2P1E <- FeaturePlot(macrophages_pf_tx_seurat, "VSIG4", split.by = "Tissue", label = T)

figR2P1AB <- ggarrange(figR2P1A, figR2P1B, nrow = 1, ncol = 2, labels = c("A", "B"))
figR2P1CDE <- ggarrange(figR2P1C, figR2P1D, figR2P1E[[1]], figR2P1E[[2]], nrow = 1, ncol = 4, labels = c("C", "D", "E", ""), widths = c(1.5, 0.5, 0.5, 0.5))

figR2P1 <- ggarrange(figR2P1AB, figR2P1CDE, nrow = 2, ncol = 1)

ggsave(filename = "figR2P1.pdf", width = 22.5, height = 10, units = "in")
