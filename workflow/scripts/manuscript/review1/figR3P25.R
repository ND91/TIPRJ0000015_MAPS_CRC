#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 25.

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(pheatmap)
require(ggpubr)

seurat_rds <- "output/q2_crc_vs_hc/subsets/hc_crcpmp_pf_SeuratObject.Rds"

group_colors <- c(`HC` = "#93CEC1", `CRC+` = "#996633")

seuratObject <- readRDS(seurat_rds)

figR3P25A <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l2, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  dplyr::filter(manual_l2 %in% c("Macrophages")) %>%
  ggplot(aes(x = manual_l2, y = Ncellperc)) +
  geom_boxplot(outlier.shape = NA, aes(col = Group)) +
  geom_point(position = position_dodge(width=0.75), aes(group = Group)) +
  labs(y = "%CD45+") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

figR3P25B <- seuratObject@meta.data %>%
  dplyr::group_by(manual_l3, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  dplyr::filter(manual_l3 %in% c("Macrophages C1Q+", "Macrophages SPP1+", "Macrophages VCAN+", "Macrophages VCAN+C1Q+")) %>%
  ggplot(aes(x = manual_l3, y = Ncellperc)) +
  geom_boxplot(outlier.shape = NA, aes(col = Group)) +
  geom_point(position = position_dodge(width=0.75), aes(group = Group)) +
  labs(y = "%CD45+") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

figR3P25 <- ggarrange(figR3P25A, figR3P25B, nrow = 1, ncol = 2, labels = c("A", "B"), widths = c(1, 4), align = "hv", legend = "bottom", common.legend = T)
ggsave(filename = "figR3P25.pdf", width = 7.5, height = 5, units = "in")
