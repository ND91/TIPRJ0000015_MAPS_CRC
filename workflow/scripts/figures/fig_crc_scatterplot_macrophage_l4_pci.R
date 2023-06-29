#!/usr/bin/env Rscript
# The goal of this script is to create a scatterplot regressing macrophage l4 (relative to all macrophages) against PCI.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_xlsx <- args[2]
scatterplot_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)

seuratObject_macrophages <- seuratObject[,seuratObject@meta.data$manual_l2 == "Macrophages" & seuratObject@meta.data$Tissue %in% c("PF", "TX")]

proportions_df <- seuratObject_macrophages@meta.data %>%
  dplyr::mutate(manual_l4 = factor(manual_l4)) %>%
  dplyr::group_by(manual_l4, SampleID, Tissue, .drop = F) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(manual_l4), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(seuratObject_macrophages@meta.data %>%
                     dplyr::select(SampleID, PCI) %>%
                     unique(),
                   by = "SampleID")

ggplot(proportions_df, aes(x = PCI, y = Ncellprop, col = Tissue)) +
  geom_smooth(method=lm) +
  geom_point() +
  theme_bw() +
  facet_wrap(~manual_l4) +
  theme(legend.pos = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 6, height = 3.25, file = dotplot_pdf)
print(dotplot)
dev.off()

# Boxplot

markerexpr_long <- data.frame(GeneID = rownames(marker_gex),
                              marker_gex) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "expr") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l4), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = unique(marker_genes$celltype)),
                GeneID = factor(GeneID, levels = markergenes_of_interest))

boxplotobj <- markerexpr_long %>%
  dplyr::filter(expr != 0) %>%
  ggplot(aes(x = Celltype, y = expr, fill = Celltype)) +
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

pdf(width = 3.5, height = 12.5, file = boxplot_pdf, bg = "white")
print(boxplotobj)
dev.off()

# Violinplot

violinplotobj_v1 <- lapply(VlnPlot(seuratObject, markergenes_of_interest, ncol = 1, combine = F), function(i){
  i + 
    theme_bw() +
    scale_fill_manual(values = macrophages_colors) +
    theme(legend.pos = "none", 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0, "lines"))
})

violinplotobj_v1 <- c(lapply(violinplotobj_v1[1:(length(violinplotobj_v1)-1)], function(i){
  i + theme(axis.text.x = element_blank())
}), violinplotobj_v1[length(violinplotobj_v1)])

pdf(width = 3.5, height = 30, file = violinplot_pdf, bg = "white")
print(ggarrange(plotlist = violinplotobj_v1, ncol = 1, ))
dev.off()

violinplotobj_v2 <- markerexpr_long %>%
  dplyr::filter(#expr != 0,
    GeneID == "CSF1R") %>%
  ggplot(aes(x = Celltype, y = expr, fill = Celltype)) +
  #geom_jitter_rast(alpha = 0.1) +
  geom_violin() +
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

pdf(width = 3.5, height = 12.5, file = violinplot_pdf, bg = "white")
print(violinplotobj_v2)
dev.off()
