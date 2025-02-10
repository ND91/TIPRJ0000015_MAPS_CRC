#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 4, point 4.

require(Seurat)
require(SingleCellExperiment)
require(TSCAN)
require(dplyr)
require(ggplot2)
require(ggrastr)

sce_hc_pf_macrophages_trajectory_rds <- "output/q1_pf_characterization/subsets/hc_pf_macrophages_trajectory_sce.Rds"
sce_crcpmp_tx_macrophages_trajectory_rds <- "output/q3_pm_tx_characterization/subsets/crcpmp_tx_macrophages_trajectory_sce.Rds"
hc_crcpmp_pf_tx_macrophages_seurat_rds <- "output/q3_pm_tx_characterization/subsets/hc_crcpmp_pf_tx_macrophages_SeuratObject.Rds"
celltype_markers_xlsx <- "config/order/celltype_markers.xlsx"

boxplot_hc_pf_macrophages_entropy_pdf <- "boxplot_hc_pf_macrophages_entropy.pdf"
boxplot_hc_crcpmp_pf_tx_macrophages_entropy_pdf <- "boxplot_hc_crcpmp_pf_tx_macrophages_entropy.pdf"

hc_crcpmp_pf_tx_macrophages_seurat <- readRDS(hc_crcpmp_pf_tx_macrophages_seurat_rds)
sce_hc_pf_macrophages_trajectory <- readRDS(sce_hc_pf_macrophages_trajectory_rds)
sce_crcpmp_tx_macrophages_trajectory <- readRDS(sce_crcpmp_tx_macrophages_trajectory_rds)

celltype_order_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(colData(sce_hc_pf_macrophages_trajectory)[,"manual_l3"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l3_colors <- celltype_order_l3$color
names(manual_l3_colors) <- celltype_order_l3$celltype

# HC-PF and CRCPMP-TX

boxplot_hc_pf_macrophages_entropy_ggplotobj <- colData(sce_hc_pf_macrophages_trajectory) %>%
  data.frame() %>%
  dplyr::select(CellID, Tissue, manual_l3, entropy) %>%
  dplyr::rows_append(colData(sce_crcpmp_tx_macrophages_trajectory) %>%
                       data.frame() %>%
                       dplyr::select(CellID, Tissue, manual_l3, entropy)) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -entropy), y = entropy, fill = manual_l3)) +
  geom_point_rast(position = position_jitterdodge(), alpha = 0.3, size = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Entropy") +
  facet_grid(.~Tissue) +
  theme_bw() +
  scale_fill_manual(values = manual_l3_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# HC-PF, CRCPMP-PF, CRCPMP-TX

sce_hc_crcpmp_pf_tx_macrophages <- as.SingleCellExperiment(hc_crcpmp_pf_tx_macrophages_seurat)
colData(sce_hc_crcpmp_pf_tx_macrophages)$entropy <- perCellEntropy(sce_hc_crcpmp_pf_tx_macrophages)

boxplot_hc_crcpmp_pf_tx_macrophages_entropy_ggplotobj <- colData(sce_hc_crcpmp_pf_tx_macrophages) %>%
  data.frame() %>%
  dplyr::mutate(Grouping = factor(paste0(Group, " ", Tissue), levels = c("HC PF", "CRC+ PF", "CRC+ TX"))) %>%
  dplyr::select(CellID, Grouping, manual_l3, entropy) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -entropy), y = entropy, fill = manual_l3)) +
  geom_point_rast(position = position_jitterdodge(), alpha = 0.3, size = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Entropy") +
  facet_grid(.~Grouping) +
  theme_bw() +
  scale_fill_manual(values = manual_l3_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(boxplot_hc_crcpmp_pf_tx_macrophages_entropy_pdf, width = 9, height = 4)
