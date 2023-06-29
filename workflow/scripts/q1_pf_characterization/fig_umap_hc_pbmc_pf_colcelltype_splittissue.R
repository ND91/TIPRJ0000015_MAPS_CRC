#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMC and PF samples from HC colored by celltype, split for tissue.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
celltype_markers_xlsx <- args[2]
celltype_level <- args[3]
hc_umap_pbmc_pf_colcelltype_splittissue_pdf <- args[4]

seuratObject <- readRDS(seurat_rds)
celltype_order <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == celltype_level,
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seuratObject@meta.data[,celltype_level])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

celltype_colors_list <- celltype_order$color
names(celltype_colors_list) <- celltype_order$number_subset

umap_df <- data.frame(CB = colnames(seuratObject),
                      Embeddings(seuratObject[["umap"]]),
                      seuratObject@meta.data) %>%
  dplyr::mutate(celltype = factor(UQ(rlang::sym(celltype_level)), levels = celltype_order$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order$number_subset))

hc_umap_pbmc_pf_colcelltype_splittissue <- umap_df %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, col = celltype_w_number)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y, fill = celltype_w_number),
                   alpha = 0.7, 
                   show.legend = F,
                   col = "#000000") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_colors_list) +
  scale_fill_manual(values = celltype_colors_list) +
  facet_wrap(~Tissue, ncol = 2) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf(width = 10, height = 6, file = hc_umap_pbmc_pf_colcelltype_splittissue_pdf)
print(hc_umap_pbmc_pf_colcelltype_splittissue)
dev.off()

sessionInfo()