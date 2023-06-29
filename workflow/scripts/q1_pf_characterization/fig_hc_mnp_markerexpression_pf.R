#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMC and PF samples from HC colored by celltype, split for tissue.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds <- args[1]
hc_violinplot_mnpmarkexpression_pdf <- args[2]
hc_boxplot_mnpmarkexpression_pdf <- args[3]
hc_dotplot_mnpmarkexpression_pdf <- args[4]
hc_umap_mnpmarkexpression_pdf <- args[5]

#celltype_gex_markers <- c("CD14", "FCGR3A", "CD1C", "CD163", "IL6", "VCAN", "TGFB1", "IL10", "C1QA", "SPP1", "MAF")
celltype_gex_markers <- c("CD1C", "CLEC9A", "XCR1", "CADM1", "VCAN", "C1QA", "CD163", "FCGR3A", "CLEC4C", "IL3RA")
celltype_pex_markers <- c("Hu.CD14", "Hu.CD16", "Hu.CD163", "Hu.CD64")
celltype_order <- c("Classical monocytes", "M1-like", "M2-like", "M1/M2-like", "CDC1s", "CDC2s", "PDCs")

seuratObject <- readRDS(seurat_rds)
seuratObject <- NormalizeData(seuratObject, assay = "CITE")

seuratObject <- seuratObject[,seuratObject@meta.data$manual_l3 != "Non-classical monocytes"]

# Violinplot

gex <- GetAssayData(seuratObject, assay = "RNA")
pex <- GetAssayData(seuratObject, assay = "CITE")

ol_genes <- celltype_gex_markers[celltype_gex_markers %in% rownames(gex)]
ol_proteins <- celltype_pex_markers[celltype_pex_markers %in% rownames(pex)]

gex_markergenes_df <- gex[ol_genes,]
pex_markergenes_df <- pex[ol_proteins,]

markerexpr_long <- data.frame(FeatureID = c(ol_genes, ol_proteins),
                              rbind(gex_markergenes_df, pex_markergenes_df)) %>%
  tidyr::pivot_longer(-FeatureID, names_to = "CB", values_to = "expr") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = celltype_order),
                FeatureID = factor(FeatureID, levels = c(ol_genes, ol_proteins)))

celltype_features_proper <- markerexpr_long %>%
  dplyr::filter(expr != 0)%>%
  dplyr::group_by(Celltype, FeatureID) %>%
  dplyr::summarize(counts = n()) %>%
  dplyr::filter(counts>2) %>%
  dplyr::mutate(celltype_featureid = paste0(Celltype, FeatureID))

violinplotobj <- markerexpr_long %>%
  dplyr::mutate(celltype_featureid = paste0(Celltype, FeatureID)) %>%
  dplyr::filter(celltype_featureid %in% celltype_features_proper$celltype_featureid) %>%
  ggplot(aes(x = Celltype, y = expr, fill = Celltype)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA) +
  facet_grid(FeatureID~., scales = "free_y") +
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

pdf(width = 3.5, height = 12.5, file = hc_violinplot_mnpmarkexpression_pdf, bg = "white")
print(violinplotobj)
dev.off()

# Boxplot

boxplotobj <- markerexpr_long %>%
  dplyr::filter(expr != 0) %>%
  ggplot(aes(x = Celltype, y = expr, fill = Celltype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(FeatureID~., scales = "free_y") +
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

pdf(width = 3.5, height = 12.5, file = hc_boxplot_mnpmarkexpression_pdf, bg = "white")
print(boxplotobj)
dev.off()

# Dotplot

markerexpr_median <- data.frame(FeatureID = c(ol_genes, ol_proteins),
                                rbind(gex_markergenes_df, pex_markergenes_df)) %>%
  tidyr::pivot_longer(-FeatureID, names_to = "CB", values_to = "exprs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(exprs)),
                   Percentage = mean(exprs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, celltype_order),
                FeatureID = factor(FeatureID, c(ol_genes, ol_proteins)))

dotplot <- markerexpr_median %>% 
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = FeatureID, y = Celltype, col = Median)) +
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

pdf(width = 7, height = 4, file = hc_dotplot_mnpmarkexpression_pdf, bg = "white")
print(dotplot)
dev.off()

# UMAP

pdf(width = 20, height = 10, file = hc_umap_mnpmarkexpression_pdf, bg = "white")
DimPlot(seuratObject, group.by = "manual_l3", label = T) + theme(legend.pos = "none") + 
  FeaturePlot(seuratObject, c(ol_genes, ol_proteins))
dev.off()

sessionInfo()