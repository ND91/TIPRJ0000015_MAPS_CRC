#!/usr/bin/env Rscript

# The goal of this script is to create a violinplot of some markers of interest.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
celltype_markers_xlsx <- args[2]
celltype_level <- args[3]
hc_violinplot_pbmc_pf_pdf <- args[4]

seuratObject <- readRDS(seurat_rds)
celltype_order <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% unique(seuratObject@meta.data[,celltype_level]),
                !is.na(canonical_marker)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::filter(row_number()==1)

counts <- GetAssayData(seuratObject, assay = "RNA")

gex_markergenes_df <- counts[rownames(counts) %in% celltype_order$canonical_marker,]

gex_markergenes_long <- data.frame(GeneID = rownames(gex_markergenes_df), gex_markergenes_df) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l2), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = unique(celltype_order$celltype)),
                GeneID = factor(GeneID, levels = unique(celltype_order$canonical_marker)))

pdf(width = 5, height = 15, file = hc_violinplot_pbmc_pf_pdf, bg = "white")
gex_markergenes_long %>%
  dplyr::filter(nUMIs != 0) %>%
  ggplot(aes(x = Celltype, y = nUMIs, fill = Celltype)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(GeneID~., scales = "free_y") +
  theme_bw() +
  theme(legend.pos = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0, "lines"))
dev.off()

sessionInfo()