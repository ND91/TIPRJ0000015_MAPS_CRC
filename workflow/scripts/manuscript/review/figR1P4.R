#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 4.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(ggpubr)
require(ggplotify)
require(ggrastr)
require(viridis)

seurat_rds <- args[1] #"output/subsets/live_singlet_nonproliferating_SeuratObject.Rds"

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Group %in% c("CRC+", "HC"),
                Tissue %in% c("PF"),
                manual_l2 %in% c("Macrophages")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

# Normalize, reduce, and recluster
seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]
options(future.globals.maxSize= 20000*1024^2)
seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 4132)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:47)
seuratObject <- FindClusters(seuratObject, resolution = 1, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:47, seed.use = 5132)

figR1P4A <- DimPlot(seuratObject) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))
figR1P4A <- LabelClusters(plot = figR1P4A, id = 'ident', box = T)

figR1P4B <- DimPlot(seuratObject, group.by = "SampleID") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))
figR1P4B <- LabelClusters(plot = figR1P4B, id = 'SampleID', box = T)

figR1P4C <- as.ggplot(pheatmap(prop.table(table(seuratObject@meta.data$SampleID, seuratObject@meta.data$seurat_clusters), margin = 2), 
                               cluster_cols = F, 
                               cluster_rows = F, 
                               display_numbers = T))

figR1P4D_list <- lapply(
  data.frame(CB = colnames(seuratObject),
             Embeddings(seuratObject[["umap"]]),
             seuratObject@meta.data,
             expr = GetAssayData(seuratObject, layers = "data")["MARCO",],
             Feature = "MARCO") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["CD163",],
                                  Feature = "CD163")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["VCAN",],
                                  Feature = "VCAN")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["CCR2",],
                                  Feature = "CCR2")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["SPP1",],
                                  Feature = "SPP1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["C1QA",],
                                  Feature = "C1QA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["CD14",],
                                  Feature = "CD14")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["FCGR3A",],
                                  Feature = "FCGR3A")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                  Feature = as.factor(Feature)) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% dplyr::arrange(expr) %>%
        ggplot(aes(x = umap_1, y = umap_2)) +
        geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
        geom_point_rast(data = . %>%
                          dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~Feature) +
        guides(colour = guide_legend(override.aes = list(size = 3),
                                     title = "Expression")) +
        scale_color_viridis() +
        #scale_colour_gradient(low = "#ffffff", high = "#FF0000") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              strip.text.x = element_text(face = "bold"))
    }
)
figR1P4D <- ggarrange(plotlist = figR1P4D_list, nrow = 2, ncol = 4)

figR1P4ABC <- ggarrange(figR1P4A, figR1P4B, figR1P4C, nrow = 1, ncol = 3, labels = c("A", "B", "C"))
figR1P4 <- ggarrange(figR1P4ABC, figR1P4D, nrow = 2, ncol = 1, heights = c(1, 1.5), labels = c("", "D"))

# ggsave(filename = "figR1P4.pdf", width = 20, height = 12, units = "in")
ggsave(filename = "figR1P4.pdf", width = 20, height = 16, units = "in")
