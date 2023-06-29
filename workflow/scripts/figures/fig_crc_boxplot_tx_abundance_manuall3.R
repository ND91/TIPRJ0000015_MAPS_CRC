#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrastr))
suppressPackageStartupMessages(require(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0("Script needs 8 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
boxplot_tx_top5_pdf <- args[2]
boxplot_tx_top10_pdf <- args[3]
boxplot_tx_all_pdf <- args[4]

seuratObject <- readRDS(seurat_rds_path)

proportions_df <- seuratObject@meta.data %>%
    dplyr::group_by(manual_l3, SampleID, Tissue, .drop = T) %>% 
    dplyr::summarise(Ncells = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(Nlrefsample = sum(Ncells),
                  Ncellprop = Ncells/Nlrefsample,
                  Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                  Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))

boxplot_tx_top5 <- proportions_df %>%
  dplyr::filter(manual_l3 %in% c("CD8 TEM", "CD4 Treg", "M2-like", "CD4 TEM", "CD4 TCM")) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to CD45+ in TX",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 3, height = 5, file = boxplot_tx_top5_pdf)
print(boxplot_tx_top5)
dev.off()

boxplot_tx_top10 <- proportions_df %>%
  dplyr::filter(manual_l3 %in% c("CD8 TEM", "CD4 Treg", "M2-like", "CD4 TEM", "CD4 TCM", "CD4 CTL", "CDC2s", "CD8 CTL", "M1/M2-like", "B intermediate lambda")) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to CD45+ in TX",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 6, height = 5, file = boxplot_tx_top10_pdf)
print(boxplot_tx_top10)
dev.off()

boxplot_tx_all <- proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(subtitle = "Proportion relative to CD45+ in TX",
       y = "Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 15, height = 5, file = boxplot_tx_all_pdf)
print(boxplot_tx_all)
dev.off()

sessionInfo()