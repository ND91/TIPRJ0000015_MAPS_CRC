#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 41.

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggrastr)
require(ggpubr)
require(ggrepel)

seurat_gcpmp_pf_rds <- "resources/gastric_cancer/maps_gc_pf_txpts.Rds"
seurat_gcpmp_pf_mnp_rds <- "resources/gastric_cancer/maps_gc_pf_mnp_txpts.Rds"
seurat_gcpmp_pf_macrophages_rds <- "resources/gastric_cancer/maps_gc_pf_macrophages_txpts.Rds"
seurat_gcpmp_tx_macrophages_rds <- "resources/gastric_cancer/maps_gc_tx_macrophages.Rds"

celltype_order_l3_monocytes_macrophages <- data.frame(celltype = c("Classical monocytes", "Non-classical monocytes", "Macrophages VCAN+", "Macrophages C1Q+", "Macrophages VCAN+C1Q+", "Macrophages SPP1+"),
                                                      color = c("#FF7F00", "#FFFF33", "#FDB462", "#B3DE69", "#BC80BD", "#FCCDE5"),
                                                      number_subset = c("1. Classical monocytes", "2. Non-classical monocytes", "3. Macrophages VCAN+", "4. Macrophages C1Q+", "5. Macrophages VCAN+C1Q+", "6. Macrophages SPP1+"))

manual_l3_monocytes_macrophages_number_colors <- celltype_order_l3_monocytes_macrophages$color
names(manual_l3_monocytes_macrophages_number_colors) <- celltype_order_l3_monocytes_macrophages$number_subset

# PF immune

seurat_gcpmp_pf_tx <- readRDS(seurat_gcpmp_pf_rds)

seurat_gcpmp_pf

# PF MNP

# PF Macrophages

seurat_gcpmp_pf_macrophages <- readRDS(seurat_gcpmp_pf_macrophages_rds)

figR3P41A <- FeaturePlot(seurat_gcpmp_pf_macrophages, features = c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), label = T, ncol = 4)

ggsave(filename = "figR3P41.pdf", width = 18, height = 8, units = "in")

# TX Macrophages

seurat_gcpmp_tx_macrophages <- readRDS(seurat_gcpmp_tx_macrophages_rds)

umap_gcpmp_tx_macrophages_ggplotobj <- data.frame(CB = colnames(seurat_gcpmp_tx_macrophages),
                                                  Embeddings(seurat_gcpmp_tx_macrophages [["umap"]]),
                                                  seurat_gcpmp_tx_macrophages@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_l3_monocytes_macrophages$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_l3_monocytes_macrophages$number_subset)) %>% 
  ggplot(aes(x = umap_1, y = umap_2)) +
  #geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 3, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = . %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = umap_1),
                               y = median(x = umap_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1,
                   size=8,
                   label.size=0.25,
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l3_monocytes_macrophages_number_colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

ggsave(filename = "figR3P41.pdf", plot = umap_gcpmp_tx_macrophages_ggplotobj, width = 5, height = 5, units = "in")
