#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(ggpubr)
require(viridis)
require(SummarizedExperiment)
require(slingshot)
require(SingleCellExperiment)

crcpmp_pf_tx_paired_macrophages_sce_ss_rds <- args[1]

crcpmp_pf_tx_paired_macrophages_sce_ss <- readRDS(crcpmp_pf_tx_paired_macrophages_sce_ss_rds)

crcpmp_pf_tx_paired_macrophages_sce_ss_sds <- SlingshotDataSet(crcpmp_pf_tx_paired_macrophages_sce_ss)
crcpmp_pf_tx_paired_macrophages_curves <- slingCurves(crcpmp_pf_tx_paired_macrophages_sce_ss_sds, as.df = TRUE)

crcpmp_pf_tx_paired_macrophages_df <- data.frame(reducedDims(crcpmp_pf_tx_paired_macrophages_sce_ss)$UMAP, 
                                                 colData(crcpmp_pf_tx_paired_macrophages_sce_ss)[,-which(colnames(colData(crcpmp_pf_tx_paired_macrophages_sce_ss)) == "slingshot")])

# Slingshot
crcpmp_pf_tx_paired_macrophages_ggplotobj <- ggplot(crcpmp_pf_tx_paired_macrophages_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = Tissue)) +
  geom_path(data = crcpmp_pf_tx_paired_macrophages_curves %>% arrange(Order), aes(group = Lineage)) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  #scale_color_manual(values = celltype_colors_pf_l3_list) +
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

crcpmp_pf_tx_paired_ggplotobj_split <- ggplot(crcpmp_pf_tx_paired_macrophages_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = manual_l3)) +
  geom_path(data = crcpmp_pf_tx_paired_macrophages_curves %>% arrange(Order), aes(group = Lineage)) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~Tissue) +
  #scale_color_manual(values = celltype_colors_pf_l3_list) +
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
