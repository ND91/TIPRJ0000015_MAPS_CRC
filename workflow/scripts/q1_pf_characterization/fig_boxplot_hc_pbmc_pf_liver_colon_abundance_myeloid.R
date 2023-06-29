#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances of manual_l3 relative to manual_l1 for our the myeloid cells found in our own PBMC and PF samples, as well as the publicly derived colon and liver samples.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
boxplot_pdf <- args[2]
boxplot_patanno_pdf <- args[3]

seuratObject <- readRDS(seurat_rds_path)

proportions_df <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = factor(manual_l3)) %>%
  dplyr::group_by(manual_l3, SampleID, Tissue, .drop = F) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))

boxplot_ggplot <- proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to myeloid",
       y = "Proportion") +
  theme_bw() +
  scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

boxplot_patanno_ggplot <- proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
  labs(subtitle = "Proportion relative to myeloid",
       y = "Proportion") +
  theme_bw() +
  scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 10, height = 5, file = boxplot_pdf)
print(boxplot_ggplot)
dev.off()

pdf(width = 10, height = 5, file = boxplot_patanno_pdf)
print(boxplot_patanno_ggplot)
dev.off()


sessionInfo()