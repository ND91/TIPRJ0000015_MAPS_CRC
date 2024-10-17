#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 3, point 26.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 arguments. Current input is:", args))
}

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(pheatmap)
require(ggpubr)

seurat_rds <- args[1] #"output/q3_pm_tx_characterization/subsets/crcpmp_tx_immune_SeuratObject.Rds"
celltype_markers_xlsx <- args[2] #"config/order/celltype_markers.xlsx"
pf_heatmap_order_xlsx <- args[3] #"config/order/q1_pf_heatmap_order.xlsx"

tx_seurat <- readRDS(seurat_rds)
Idents(tx_seurat) <- "manual_l3"

celltype_order_tx_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(tx_seurat@meta.data[,"manual_l3"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l3_colors_tx <- celltype_order_tx_l3$color
names(manual_l3_colors_tx) <- celltype_order_tx_l3$celltype


## Version 1

figR3P34A <- DimPlot(tx_seurat, group.by = "manual_l3", label = T) + labs(title = "PM")# + theme(legend.position = "bottom")
figR3P34B_list <- FeaturePlot(tx_seurat, c("MARCO", "CD163", "VCAN", "CCR2", "SPP1", "C1QA", "CD14", "FCGR3A"), label = F, combine = FALSE)
figR3P34B <- ggarrange(plotlist = figR3P34B_list, nrow = 2, ncol = 4)

## Version 2

figR3P34B_list <- lapply(
  data.frame(CB = colnames(tx_seurat),
             Embeddings(tx_seurat[["umap"]]),
             tx_seurat@meta.data,
             expr = GetAssayData(tx_seurat, layers = "data")["MARCO",],
             Feature = "MARCO") %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["CD163",],
                                  Feature = "CD163")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["VCAN",],
                                  Feature = "VCAN")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["CCR2",],
                                  Feature = "CCR2")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["SPP1",],
                                  Feature = "SPP1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["C1QA",],
                                  Feature = "C1QA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["CD14",],
                                  Feature = "CD14")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["FCGR3A",],
                                  Feature = "FCGR3A")) %>%
    dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                  Embeddings(tx_seurat[["umap"]]),
                                  tx_seurat@meta.data,
                                  expr = GetAssayData(tx_seurat, layers = "data")["VSIG4",],
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
figR3P34B <- ggarrange(plotlist = figR3P34B_list, nrow = 2, ncol = 5)

marker_gex <- GetAssayData(tx_seurat, assay = "RNA")[which(rownames(GetAssayData(tx_seurat, assay = "RNA")) %in% c("CCR2", "CD14", "FCGR3A", "VSIG4", "MARCO", "CD163", "VCAN", "C1QA", "SPP1")), ]

figR3P34Cv1 <- data.frame(FeatureID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-c(FeatureID), names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(tx_seurat@meta.data), 
                               Celltype = tx_seurat@meta.data[,"manual_l3"]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, celltype_order_tx_l3$celltype),
                FeatureID = factor(FeatureID, rev(c("CCR2", "CD14", "FCGR3A", "VSIG4", "MARCO", "CD163", "VCAN", "C1QA", "SPP1")))) %>% 
  ggplot(aes(x = Celltype, y = FeatureID, col = Median)) +
  geom_tile(alpha = 0, col = "grey") +
  geom_point(aes(size = Percentage)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x=unit(0, "lines"))

print(figR3P34Cv1)
ggsave(filename = "figR3P34Cv1.pdf", width = 7.5, height = 4.25, units = "in")

figR3P34Cv2 <- data.frame(FeatureID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-c(FeatureID), names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(tx_seurat@meta.data), 
                               Celltype = tx_seurat@meta.data[,"manual_l3"]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::filter(Celltype %in% c("Classical monocytes", "Non-classical monocytes", "Macrophages VCAN+", "Macrophages C1Q+", "Macrophages VCAN+C1Q+", "Macrophages SPP1+")) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, c("Classical monocytes", "Non-classical monocytes", "Macrophages VCAN+", "Macrophages C1Q+", "Macrophages VCAN+C1Q+", "Macrophages SPP1+")),
                FeatureID = factor(FeatureID, rev(c("CCR2", "CD14", "FCGR3A", "VSIG4", "MARCO", "CD163", "VCAN", "C1QA", "SPP1")))) %>% 
  ggplot(aes(x = Celltype, y = FeatureID, col = Median)) +
  geom_tile(alpha = 0, col = "grey") +
  geom_point(aes(size = Percentage)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x=unit(0, "lines"))

print(figR3P34Cv2)
ggsave(filename = "figR3P34Cv2.pdf", width = 2, height = 4.25, units = "in")

figR3P34Dv1 <- data.frame(CB = colnames(tx_seurat),
                          Embeddings(tx_seurat[["umap"]]),
                          tx_seurat@meta.data,
                          expr = GetAssayData(tx_seurat, assay = "RNA")["CCR2",],
                          Feature = "CCR2") %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["CD14",],
                                Feature = "CD14")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["FCGR3A",],
                                Feature = "FCGR3A")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["VSIG4",],
                                Feature = "VSIG4")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["MARCO",],
                                Feature = "MARCO")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["CD163",],
                                Feature = "CD163")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["VCAN",],
                                Feature = "VCAN")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["C1QA",],
                                Feature = "C1QA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["SPP1",],
                                Feature = "SPP1")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                Feature = factor(Feature, levels = c("CCR2", "CD14", "FCGR3A", "VSIG4", "MARCO", "CD163", "VCAN", "C1QA", "SPP1")),
                celltype = factor(manual_l3, levels = celltype_order_tx_l3$celltype)) %>%
  ggplot(aes(x = celltype, y = expr, fill = celltype)) +
  geom_violin(trim = F, scale = "width") +
  facet_wrap(~Feature, ncol = 1, scales = "free_y", strip.position = "right") +
  labs(y = "nUMIs") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_fill_manual(values = manual_l3_colors_tx) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "lines"))

ggsave(filename = "figR3P34Dv1.pdf", width = 7.5, height = 7.5, units = "in")

figR3P34Dv2 <- data.frame(CB = colnames(tx_seurat),
                          Embeddings(tx_seurat[["umap"]]),
                          tx_seurat@meta.data,
                          expr = GetAssayData(tx_seurat, assay = "RNA")["CCR2",],
                          Feature = "CCR2") %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["CD14",],
                                Feature = "CD14")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["FCGR3A",],
                                Feature = "FCGR3A")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["VSIG4",],
                                Feature = "VSIG4")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["MARCO",],
                                Feature = "MARCO")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["CD163",],
                                Feature = "CD163")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["VCAN",],
                                Feature = "VCAN")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["C1QA",],
                                Feature = "C1QA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(tx_seurat),
                                Embeddings(tx_seurat[["umap"]]),
                                tx_seurat@meta.data,
                                expr = GetAssayData(tx_seurat, assay = "RNA")["SPP1",],
                                Feature = "SPP1")) %>%
  dplyr::filter(manual_l3 %in% c("Classical monocytes", "Non-classical monocytes", "Macrophages VCAN+", "Macrophages C1Q+", "Macrophages VCAN+C1Q+", "Macrophages SPP1+")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                Feature = factor(Feature, levels = c("CCR2", "CD14", "FCGR3A", "VSIG4", "MARCO", "CD163", "VCAN", "C1QA", "SPP1")),
                celltype = factor(manual_l3, levels = c("Classical monocytes", "Non-classical monocytes", "Macrophages VCAN+", "Macrophages C1Q+", "Macrophages VCAN+C1Q+", "Macrophages SPP1+"))) %>%
  ggplot(aes(x = celltype, y = expr, fill = celltype)) +
  geom_violin(trim = F, scale = "width") +
  facet_wrap(~Feature, ncol = 1, scales = "free_y", strip.position = "right") +
  labs(y = "nUMIs") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_fill_manual(values = manual_l3_colors_tx) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "lines"))

ggsave(filename = "figR3P34Dv2.pdf", width = 3, height = 7.5, units = "in")

figR3P34 <- ggarrange(figR3P34A, figR3P34B, nrow = 1, ncol = 2, labels = c("A", "B"), widths = c(1.5, 2))
ggsave(filename = "figR3P34.pdf", width = 28, height = 8, units = "in")
