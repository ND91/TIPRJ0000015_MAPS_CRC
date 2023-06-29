#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of all macrophage subsets (l4) relative to the total macrophages in PF.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
boxplot_pf_pdf <- args[2]

seuratObject <- readRDS(seurat_rds_path)

proportions_df <- seuratObject@meta.data %>%
  dplyr::filter(manual_l2 == "Macrophages",
                Tissue == "PF") %>%
  dplyr::group_by(manual_l4, SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(manual_l4), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(unique(seuratObject@meta.data[,c("SampleID", "Group")]), by = "SampleID")


boxplot_pf <- proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l4, -Ncellprop), y = Ncellprop, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  #facet_grid(.~, scales = "free", space='free') +
  labs(subtitle = "Proportion relative to all macrophages in PF",
       y = "Proportion") +
  scale_color_manual(values = c(`HC` = "#91d1c2", `CRC+` = "#845422")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 5, height = 5, file = boxplot_pf_pdf)
print(boxplot_pf)
dev.off()

sessionInfo()