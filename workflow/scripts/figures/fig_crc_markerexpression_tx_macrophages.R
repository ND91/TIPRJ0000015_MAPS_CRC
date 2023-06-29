#!/usr/bin/env Rscript
# The goal of this script is to create a boxplot of l1 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_xlsx <- args[2]
markergene_analysis_rds <- args[3]
dotplot_pdf <- args[4]
boxplot_pdf <- args[5]
violinplot_pdf <- args[6]

seuratObject <- readRDS(seurat_rds)

markergene_analysis <- readRDS(markergene_analysis_rds)

macrophages_colors_df <- unique(marker_genes[,c("celltype", "color")])
macrophages_colors <- macrophages_colors_df$color
names(macrophages_colors) <- macrophages_colors_df$celltype

marker_genes <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == "manual_l4",
                celltype %in% unique(seuratObject@meta.data[,"manual_l4"]),
                modality == "gene")

markergenes_of_interest <- c("CSF1R", "MARCO", "MSR1", "MRC1", "CD163", "CXCR1", "CXCR4", "HLA-DRA", "IL10", "TGFB1", "ARG1", "VEGFA", "MAF", "RETNLB")
markergenes_of_interest <- markergenes_of_interest[markergenes_of_interest %in% rownames(seuratObject)]

marker_gex <- GetAssayData(seuratObject, assay = "RNA")[which(rownames(GetAssayData(seuratObject, assay = "RNA")) %in% markergenes_of_interest), ]

marker_gex_median <- data.frame(GeneID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data[,"manual_l4"]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, unique(marker_genes$celltype)),
                GeneID = factor(GeneID, markergenes_of_interest))

dotplot <- marker_gex_median %>% 
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = GeneID, y = Celltype, col = Median)) +
  geom_tile(alpha = 0, col = "grey") +
  geom_point(aes(size = Percentage)) +
  theme_bw() +
  theme(legend.pos = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x=unit(0, "lines"))

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
