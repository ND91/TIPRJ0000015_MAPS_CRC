#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 2, point 6.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggrastr)
require(viridis)
require(pheatmap)
require(ggpubr)

seurat_rds <- args[1] #"output/q1_pf_characterization/subsets/hc_pf_SeuratObject.Rds"
seurat_myeloid_rds <- args[2] #"output/q1_pf_characterization/subsets/hc_pf_Myeloid_SeuratObject.Rds"

seuratObject_pf <- readRDS(seurat_rds)
seuratObject_pf_myeloid <- readRDS(seurat_myeloid_rds)

figR2P6A_list <- lapply(
  data.frame(CB = colnames(seuratObject_pf),
             Embeddings(seuratObject_pf[["wnn.umap"]]),
             seuratObject_pf@meta.data,
             expr = GetAssayData(seuratObject_pf, layers = "data", assay = "CITE")["Hu.CD194",],
             Feature = "CD194") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf),
                                  Embeddings(seuratObject_pf[["wnn.umap"]]),
                                  seuratObject_pf@meta.data,
                                  expr = GetAssayData(seuratObject_pf, layers = "data", assay = "CITE")["Hu.CD195",],
                                  Feature = "CD195")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf),
                                  Embeddings(seuratObject_pf[["wnn.umap"]]),
                                  seuratObject_pf@meta.data,
                                  expr = GetAssayData(seuratObject_pf, layers = "data", assay = "CITE")["Hu.CD196",],
                                  Feature = "CD196")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                  Feature = as.factor(Feature)) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% dplyr::arrange(expr) %>%
        ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
        geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
        geom_point_rast(data = . %>%
                          dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~Feature) +
        guides(colour = guide_legend(override.aes = list(size = 3),
                                     title = "Expression")) +
        scale_colour_gradient(low = "#ffffff", high = "#FF0000") +
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

figR2P6B_list <- lapply(
  data.frame(CB = colnames(seuratObject_pf),
             Embeddings(seuratObject_pf[["wnn.umap"]]),
             seuratObject_pf@meta.data,
             expr = GetAssayData(seuratObject_pf, layers = "data", assay = "RNA")["IL12A",],
             Feature = "IL12A") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf),
                                  Embeddings(seuratObject_pf[["wnn.umap"]]),
                                  seuratObject_pf@meta.data,
                                  expr = GetAssayData(seuratObject_pf, layers = "data", assay = "RNA")["IL1A",],
                                  Feature = "IL1A")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                  Feature = as.factor(Feature)) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% dplyr::arrange(expr) %>%
        ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
        geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
        geom_point_rast(data = . %>%
                          dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~Feature) +
        guides(colour = guide_legend(override.aes = list(size = 3),
                                     title = "Expression")) +
        # scale_colour_gradient(low = "#ffffff", high = "#FF0000") +
        scale_color_viridis() +
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

figR2P6C_list <- lapply(
  data.frame(CB = colnames(seuratObject_pf_myeloid),
             Embeddings(seuratObject_pf_myeloid[["wnn.umap"]]),
             seuratObject_pf_myeloid@meta.data,
             expr = GetAssayData(seuratObject_pf_myeloid, assay = "CITE")["Hu.CD194",],
             Feature = "CD194") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf_myeloid),
                                  Embeddings(seuratObject_pf_myeloid[["wnn.umap"]]),
                                  seuratObject_pf_myeloid@meta.data,
                                  expr = GetAssayData(seuratObject_pf_myeloid, assay = "CITE")["Hu.CD195",],
                                  Feature = "CD195")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf_myeloid),
                                  Embeddings(seuratObject_pf_myeloid[["wnn.umap"]]),
                                  seuratObject_pf_myeloid@meta.data,
                                  expr = GetAssayData(seuratObject_pf_myeloid, assay = "CITE")["Hu.CD196",],
                                  Feature = "CD196")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first")) %>%
    dplyr::filter(manual_l2 %in% c("Macrophages", "CDCs", "PDCs")) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% 
        dplyr::arrange(expr) %>%
        ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
        geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
        geom_point_rast(data = . %>%
                          dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~Feature) +
        guides(colour = guide_legend(override.aes = list(size = 3),
                                     title = "Expression")) +
        xlim(-6,11) +
        ylim(-5,8) +
        scale_colour_gradient(low = "#ffffff", high = "#FF0000") +
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

figR2P6D_list <- lapply(
  data.frame(CB = colnames(seuratObject_pf_myeloid),
             Embeddings(seuratObject_pf_myeloid[["wnn.umap"]]),
             seuratObject_pf_myeloid@meta.data,
             expr = GetAssayData(seuratObject_pf_myeloid, layers = "data", assay = "RNA")["IL12A",],
             Feature = "IL12A") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject_pf_myeloid),
                                  Embeddings(seuratObject_pf_myeloid[["wnn.umap"]]),
                                  seuratObject_pf_myeloid@meta.data,
                                  expr = GetAssayData(seuratObject_pf_myeloid, layers = "data", assay = "RNA")["IL1A",],
                                  Feature = "IL1A")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                  Feature = as.factor(Feature)) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% dplyr::arrange(expr) %>%
        ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
        geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
        geom_point_rast(data = . %>%
                          dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~Feature) +
        guides(colour = guide_legend(override.aes = list(size = 3),
                                     title = "Expression")) +
        xlim(-6,11) +
        ylim(-5,8) +
        # scale_colour_gradient(low = "#ffffff", high = "#FF0000") +
        scale_color_viridis() +
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

figR2P6A <- ggarrange(plotlist = figR2P6A_list, nrow = 1, ncol = 3, legend = "bottom", common.legend = T)
figR2P6B <- ggarrange(plotlist = figR2P6B_list, nrow = 1, ncol = 3, legend = "bottom", common.legend = T)
figR2P6C <- ggarrange(plotlist = figR2P6C_list, nrow = 1, ncol = 3, legend = "bottom", common.legend = T)
figR2P6D <- ggarrange(plotlist = figR2P6D_list, nrow = 1, ncol = 3, legend = "bottom", common.legend = T)

figR2P6 <- ggarrange(figR2P6A, figR2P6B, figR2P6C, figR2P6D, nrow = 4, labels = c("A) Immune CITE", "B) Immune GEX", "C) MNP CITE", "D) MNP GEX"))

# ggsave(filename = "figR2P6.pdf", width = 18, height = 8, units = "in")
ggsave(filename = "figR2P6.pdf", width = 9, height = 12, units = "in")

