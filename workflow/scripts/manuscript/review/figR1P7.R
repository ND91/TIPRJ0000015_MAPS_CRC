#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 7.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggrastr)
require(ggpubr)
require(viridis)

seurat_crcpmp_tx_monocytes_macrophages_rds <- args[1] #"output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages_SeuratObject.Rds"

seuratObject <- readRDS(seurat_crcpmp_tx_monocytes_macrophages_rds)

figR1P7_list <- lapply(
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
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["VSIG4",],
                                  Feature = "VSIG4")) %>%
    dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                  Feature = as.factor(Feature)) %>%
    split(., .$Feature), FUN = function(plotdf){
      plotdf %>% dplyr::arrange(expr) %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
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
figR1P7 <- ggarrange(plotlist = figR1P7_list, nrow = 2, ncol = 5, legend = "bottom", common.legend = T)

# ggsave(filename = "figR1P7.pdf", width = 18, height = 8, units = "in")
ggsave(filename = "figR1P7.pdf", width = 20, height = 8.5, units = "in")
