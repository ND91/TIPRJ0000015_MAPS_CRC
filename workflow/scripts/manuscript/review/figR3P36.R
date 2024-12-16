#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 36.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggpubr)

seurat_rds <- args[1] #"output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds"
celltype_markers_xlsx <- args[2] #"config/order/celltype_markers.xlsx"

hc_pf_seurat <- readRDS(seurat_rds)
Idents(hc_pf_seurat) <- "manual_l3"

celltype_order_hc_pf_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(hc_pf_seurat@meta.data[,"manual_l3"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l3_colors_hc_pf <- celltype_order_hc_pf_l3$color
names(manual_l3_colors_hc_pf) <- celltype_order_hc_pf_l3$celltype

#V1
figR3P36 <- VlnPlot(hc_pf_seurat, features = "Hu.CD163", cols = manual_l3_colors_hc_pf, sort = names(manual_l3_colors_hc_pf))

ggsave(filename = "figR3P36_seurat.pdf", width = 12.5, height = 5, units = "in")

#V2
violinplot_hc_pf_cd163_ggplotobj <- 
  data.frame(CB = colnames(hc_pf_seurat),
             Embeddings(hc_pf_seurat[["wnn.umap"]]),
             hc_pf_seurat@meta.data,
             expr = GetAssayData(hc_pf_seurat, assay = "CITE")["Hu.CD163",],
             Feature = "CD163") %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l3, levels = celltype_order_hc_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_hc_pf_l3$number_subset)) %>%
  ggplot(aes(x = celltype, y = expr, fill = celltype)) +
  geom_hline(yintercept = 0) +
  geom_violin(trim = F, scale = "width") +
  geom_jitter_rast(size = 0.0001) +
  labs(title = "CD163",
       subtitle = "CITE",
       y = "nUMIs") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_fill_manual(values = manual_l3_colors_hc_pf) +
  ylim(0, NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "lines"))

print(violinplot_hc_pf_cd163_ggplotobj)

ggsave(filename = "figR3P36.pdf", plot = violinplot_hc_pf_cd163_ggplotobj, width = 9, height = 4, units = "in")

