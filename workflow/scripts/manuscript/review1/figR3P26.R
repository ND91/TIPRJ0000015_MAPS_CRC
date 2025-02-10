#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 26.

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(pheatmap)
require(ggpubr)

seurat_rds <- "output/q3_pm_tx_characterization/subsets/crcpmp_tx_monocytes_macrophages_SeuratObject.Rds"

tx_mnp_seurat <- readRDS(seurat_rds)
Idents(tx_mnp_seurat) <- "manual_l3"

## Version 1

figR3P26A <- DimPlot(tx_mnp_seurat, group.by = "manual_l3", label = T) + labs(title = "PM")# + theme(legend.position = "bottom")
figR3P26B_list <- FeaturePlot(tx_mnp_seurat, c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), label = F, combine = FALSE)
figR3P26B <- ggarrange(plotlist = figR3P26B_list, nrow = 2, ncol = 4)

## Version 2

figR3P26B_list <- lapply(
  data.frame(CB = colnames(tx_mnp_seurat),
             Embeddings(tx_mnp_seurat[["umap"]]),
             tx_mnp_seurat@meta.data,
             expr = GetAssayData(tx_mnp_seurat, layers = "data")["MARCO",],
             Feature = "MARCO") %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["CD163",],
                                  Feature = "CD163")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["VCAN",],
                                  Feature = "VCAN")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["CCR2",],
                                  Feature = "CCR2")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["SPP1",],
                                  Feature = "SPP1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["C1QA",],
                                  Feature = "C1QA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["CD14",],
                                  Feature = "CD14")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["FCGR3A",],
                                  Feature = "FCGR3A")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_mnp_seurat),
                                  Embeddings(tx_mnp_seurat[["umap"]]),
                                  tx_mnp_seurat@meta.data,
                                  expr = GetAssayData(tx_mnp_seurat, layers = "data")["VSIG4",],
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
figR3P26B <- ggarrange(plotlist = figR3P26B_list, nrow = 2, ncol = 5)

figR3P26 <- ggarrange(figR3P26A, figR3P26B, nrow = 1, ncol = 2, labels = c("A", "B"), widths = c(1.2, 2))
ggsave(filename = "figR3P26.pdf", width = 27.5, height = 8, units = "in")


