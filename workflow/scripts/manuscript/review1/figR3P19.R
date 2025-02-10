#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 19.

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggrastr)
require(pheatmap)
require(pheatmap)
require(ggpubr)
require(viridis)

seurat_rds <- "output/subsets/live_singlet_nonproliferating_SeuratObject.Rds"

seuratObject <- readRDS(seurat_rds)

selected_cells <- seuratObject@meta.data %>%
  dplyr::filter(# manual_l2 %in% c("Monocytes", "Macrophages", "CDCs"),
                manual_l3 %in% c("Classical monocytes", "Non-classical monocytes", "Macrophages C1Q+", "Macrophages SPP1+", "Macrophages VCAN+", "Macrophages VCAN+C1Q+", "CDC2s")
  )

pbmc_pf_mnp_seurat <- seuratObject[,which(seuratObject@meta.data$CellID %in% selected_cells$CellID)]

# Normalize, reduce, and recluster
pbmc_pf_mnp_seurat <- DietSeurat(pbmc_pf_mnp_seurat, counts = T, data = T, scale.data = F)
pbmc_pf_mnp_seurat <- pbmc_pf_mnp_seurat[Matrix::rowSums(pbmc_pf_mnp_seurat) != 0, ]
options(future.globals.maxSize= 20000*1024^2)
pbmc_pf_mnp_seurat <- SCTransform(pbmc_pf_mnp_seurat, conserve.memory = T)
pbmc_pf_mnp_seurat <- RunPCA(object = pbmc_pf_mnp_seurat, npcs = 100, seed.use = 3216326)
pbmc_pf_mnp_seurat <- FindNeighbors(pbmc_pf_mnp_seurat, reduction = "pca", dims = 1:48)
pbmc_pf_mnp_seurat <- FindClusters(pbmc_pf_mnp_seurat, resolution = 1, verbose = FALSE)
pbmc_pf_mnp_seurat <- RunUMAP(pbmc_pf_mnp_seurat, dims = 1:48, seed.use = 63632)

Idents(pbmc_pf_mnp_seurat) <- "manual_l3"

# Version 1

figR3P19A <- DimPlot(pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PBMC"], group.by = "manual_l3", label = T) + labs(title = "PBMC")# + theme(legend.position = "bottom")
figR3P19B_list <- FeaturePlot(pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PBMC"], c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), split.by = "Tissue", label = F, combine = FALSE)
figR3P19B <- ggarrange(plotlist = figR3P19B_list, nrow = 2, ncol = 4, align = "hv")

figR3P19AB <- ggarrange(figR3P19A, figR3P19B, nrow = 1, ncol = 2, labels = c("A", "B"), widths = c(1.2, 2))

figR3P19C <- DimPlot(pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PF"], group.by = "manual_l3", label = T) + labs(title = "PF")# + theme(legend.position = "bottom")
figR3P19D_list <- FeaturePlot(pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PF"], c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), split.by = "Tissue", label = F, combine = FALSE)
figR3P19D <- ggarrange(plotlist = figR3P19D_list, nrow = 2, ncol = 4, align = "hv")

figR3P19CD <- ggarrange(figR3P19C, figR3P19D, nrow = 1, ncol = 2, labels = c("C", "D"), widths = c(1.2, 2))

figR3P19 <- ggarrange(figR3P19AB, figR3P19CD, nrow = 2, ncol = 1)
# ggsave(filename = "figR3P19.pdf", width = 30, height = 16, units = "in")
ggsave(filename = "figR3P19.pdf", width = 25, height = 15, units = "in")

# Version 2

pbmc_mnp_seurat <- pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PBMC"]
pf_mnp_seurat <- pbmc_pf_mnp_seurat[,pbmc_pf_mnp_seurat@meta.data$Tissue == "PF"]

figR3P19A <- DimPlot(pbmc_mnp_seurat, group.by = "manual_l3", label = T) + labs(title = "PBMC")# + theme(legend.position = "bottom")
# ggsave("figR3P19Av2.pdf", plot = DimPlot(pbmc_mnp_seurat, group.by = "manual_l3", label = T, split.by = "Group") + labs(title = "PBMC"), width = 15, height = 5)
figR3P19B_list <- lapply(
  data.frame(CB = colnames(pbmc_mnp_seurat),
             Embeddings(pbmc_mnp_seurat[["umap"]]),
             pbmc_mnp_seurat@meta.data,
             expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["MARCO",],
             Feature = "MARCO") %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["CD163",],
                                  Feature = "CD163")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["VCAN",],
                                  Feature = "VCAN")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["CCR2",],
                                  Feature = "CCR2")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["SPP1",],
                                  Feature = "SPP1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["C1QA",],
                                  Feature = "C1QA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["CD14",],
                                  Feature = "CD14")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pbmc_mnp_seurat),
                                  Embeddings(pbmc_mnp_seurat[["umap"]]),
                                  pbmc_mnp_seurat@meta.data,
                                  expr = GetAssayData(pbmc_mnp_seurat, layers = "data")["FCGR3A",],
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
figR3P19B <- ggarrange(plotlist = figR3P19B_list, nrow = 2, ncol = 4)

figR3P19C <- DimPlot(pf_mnp_seurat, group.by = "manual_l3", label = T) + labs(title = "PF")# + theme(legend.position = "bottom")
figR3P19D_list <- lapply(
  data.frame(CB = colnames(pf_mnp_seurat),
             Embeddings(pf_mnp_seurat[["umap"]]),
             pf_mnp_seurat@meta.data,
             expr = GetAssayData(pf_mnp_seurat, layers = "data")["MARCO",],
             Feature = "MARCO") %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["CD163",],
                                  Feature = "CD163")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["VCAN",],
                                  Feature = "VCAN")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["CCR2",],
                                  Feature = "CCR2")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["SPP1",],
                                  Feature = "SPP1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["C1QA",],
                                  Feature = "C1QA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["CD14",],
                                  Feature = "CD14")) %>%
    dplyr::rows_append(data.frame(CB = colnames(pf_mnp_seurat),
                                  Embeddings(pf_mnp_seurat[["umap"]]),
                                  pf_mnp_seurat@meta.data,
                                  expr = GetAssayData(pf_mnp_seurat, layers = "data")["FCGR3A",],
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
figR3P19D <- ggarrange(plotlist = figR3P19D_list, nrow = 2, ncol = 4)

figR3P19AB <- ggarrange(figR3P19A, figR3P19B, nrow = 1, ncol = 2, labels = c("C", "D"), widths = c(1.2, 2))
figR3P19CD <- ggarrange(figR3P19C, figR3P19D, nrow = 1, ncol = 2, labels = c("C", "D"), widths = c(1.2, 2))

figR3P19 <- ggarrange(figR3P19AB, figR3P19CD, nrow = 2, ncol = 1)
# ggsave(filename = "figR3P19.pdf", width = 30, height = 16, units = "in")
ggsave(filename = "figR3P19.pdf", width = 25, height = 15, units = "in")
