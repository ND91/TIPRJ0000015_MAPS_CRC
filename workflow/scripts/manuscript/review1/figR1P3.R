#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 3.

BiocManager::install("UCell")

require(Seurat)
require(dplyr)
require(UCell)
require(ggplot2)
require(ggrastr)

hc_pf_T_seuratobject_rds <- "output/q1_pf_characterization/hc_pf_T_SeuratObject.Rds"
t_markers_xlsx <- "config/order/t_genesets.xlsx"

seuratObject <- readRDS(hc_pf_T_seuratobject_rds)

t_markers <- readxl::read_excel(t_markers_xlsx)
t_markers_list <- strsplit(t_markers$Markers, ";")
names(t_markers_list) <- t_markers$Celltype

DefaultAssay(seuratObject) <- "RNA"
seuratObject_annotated <- AddModuleScore_UCell(seuratObject, 
                                               features = t_markers_list, 
                                               name = NULL,
                                               ncores = 4)

seuratObject_annotated@meta.data %>%
  dplyr::select(c("CellID", "SampleID", "Donor", "manual_l3", "CD4 naive", "CD4 TCM", "CD4 Treg", "CD4 CTL", "CD4 TEM", "CD4 MAIT", "CD8 naive", "CD8 TCM", "CD8 TEM", "CD8 NKT", "CD8 ITGA1+", "CD8 MAIT", "GDT", "DNT")) %>%
  tidyr::pivot_longer(-c("CellID", "SampleID", "Donor", "manual_l3"), names_to = "Tcell", values_to = "UCell") %>%
  dplyr::group_by(manual_l3, Tcell) %>%
  # dplyr::summarize(median_UCell = median(UCell)) %>%
  dplyr::mutate(Tcell = factor(Tcell, levels = c("CD4 naive", "CD4 TCM", "CD4 Treg", "CD4 CTL", "CD4 TEM", "CD4 MAIT", "CD8 naive", "CD8 TCM", "CD8 TEM", "CD8 NKT", "CD8 ITGA1+", "CD8 MAIT", "GDT", "DNT")),
                manual_l3 = factor(manual_l3, levels = c("CD4 naive", "CD4 TCM", "CD4 Treg", "CD4 CTL", "CD4 TEM", "CD4 MAIT", "CD8 naive", "CD8 TCM", "CD8 TEM", "CD8 NKT", "CD8 ITGA1+", "CD8 MAIT", "GDT", "DNT"))) %>%
  ggplot(aes(x = manual_l3, y = UCell)) +
  # geom_jitter_rast() +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  labs(y = "UCell Score") +
  facet_grid(Tcell~.) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = "figR1P3.pdf", width = 6, height = 15, units = "in")
