#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 6.

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggrastr)
require(viridis)
require(ggpubr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1] #"output/curated/curated_SeuratObject.Rds"
celltype_markers_xlsx <- args[2] #"config/order/celltype_markers.xlsx"

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Tissue %in% "PF",
                Group %in% c("CRC+", "HC"),
                !manual_l2 %in% c("Lineage negative HP+PRG4+", "Dead/debris", "Immune SMIM25+", "Lineage negative HP+PRG4+", "Multiplets", "NK/ILC proliferating", "T apoptotic", "T proliferating", "DC proliferating", "Endothelial", "HSPCs")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

celltype_order_pf_l2_monomacs <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seuratObject@meta.data[,"manual_l2"])) %>%
  dplyr::mutate(celltype = ifelse(celltype %in% c("Monocytes", "Macrophages"), "Mono-macs", celltype),
                color = ifelse(celltype %in% c("Mono-macs"), "#E31A1C", color)) %>%
  unique() %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype)) %>%
  dplyr::mutate(color = ifelse(celltype == "HSPCs", "#D2691E", color),
                color = ifelse(celltype == "Erythroblasts", "#8B008B", color),
                color = ifelse(celltype == "Epithelial", "#4682B4", color),
                color = ifelse(celltype == "Endothelial", "#20B2AA", color)) %>%
  dplyr::rows_append(data.frame(celltype = c("Fibroblasts", "Myofibroblasts"), 
                                color = c("#9ACD32", "#8A2BE2"), 
                                number_subset = c("20. Fibroblasts", "21. Fibroblasts")))

manual_l2_monomacs_colors_pf <- celltype_order_pf_l2_monomacs$color
names(manual_l2_monomacs_colors_pf) <- celltype_order_pf_l2_monomacs$celltype


celltype_order_pf_l1l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l1",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seuratObject@meta.data[,"manual_l1"]),
                !celltype %in% c("Mesenchymal", "Myeloid")) %>%
  dplyr::rows_append(data.frame(celltype = "Mono-macs", 
                                color = "#E31A1C")) %>%
  dplyr::rows_append(data.frame(celltype = "DCs", 
                                color = "#FB9A99")) %>% 
  dplyr::rows_append(data.frame(celltype = "Granulocytes", 
                                color = "#33A02C")) %>%
  dplyr::arrange(factor(celltype, levels = c("T", "NK/ILC", "B", "Mono-macs", "DCs", "Granulocytes", "Erythroblasts", "Epithelial", "(Myo)Fibroblasts"))) %>%
  # dplyr::rows_append(data.frame(celltype = "DCs", 
  #                               color = "#FDBF6F")) %>%
  unique() %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype)) %>%
  dplyr::rows_append(data.frame(celltype = c("(Myo)Fibroblasts"), 
                                color = c("#9ACD32"), 
                                number_subset = c("10. (Myo)Fibroblasts")))

manual_l1l2_colors_pf <- celltype_order_pf_l1l2$color
names(manual_l1l2_colors_pf) <- celltype_order_pf_l1l2$celltype


seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]
options(future.globals.maxSize= 30000*1024^2)
seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 52345)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:41)
seuratObject <- FindClusters(seuratObject, resolution = 1, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:41, seed.use = 32)

seuratObject@meta.data <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l2_monomac = ifelse(manual_l2 %in% c("Monocytes", "Macrophages"), "Mono-macs", manual_l2),
                manual_l1l2 = manual_l1,
                manual_l1l2 = ifelse(manual_l2 %in% c("Myofibroblasts", "Fibroblasts"), "(Myo)Fibroblasts", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("Monocytes", "Macrophages"), "Mono-macs", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("CDCs", "PDCs"), "DCs", manual_l1l2),
                manual_l1l2 = ifelse(manual_l2 %in% c("Granulocytes"), "Granulocytes", manual_l1l2))

DimPlot(seuratObject, group.by = "manual_l3", label = T, cols = manual_l2_monomacs_colors_pf)

figR1P6A <- data.frame(CB = colnames(seuratObject),
                       Embeddings(seuratObject[["umap"]]),
                       seuratObject@meta.data) %>%
  #dplyr::mutate(celltype = manual_l2_monomac,
  dplyr::mutate(celltype = manual_l1l2,
                # celltype = factor(celltype, levels = celltype_order_pf_l2_monomacs$celltype),
                celltype = factor(celltype, levels = celltype_order_pf_l1l2$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l1l2$number_subset)) %>% 
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = . %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = umap_1),
                               y = median(x = umap_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   size = 9,
                   alpha = 1,
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  # scale_color_manual(values = manual_l2_monomacs_colors_pf) +
  scale_color_manual(values = manual_l1l2_colors_pf) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

figR1P6B_list <- lapply(
  data.frame(CB = colnames(seuratObject),
             Embeddings(seuratObject[["umap"]]),
             seuratObject@meta.data,
             expr = GetAssayData(seuratObject, layers = "data")["TGFB1",],
             Feature = "TGFB1") %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["VEGFA",],
                                  Feature = "VEGFA")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["IL10",],
                                  Feature = "IL10")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["TNF",],
                                  Feature = "TNF")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["IFNG",],
                                  Feature = "IFNG")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["IL1B",],
                                  Feature = "IL1B")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["JAK1",],
                                  Feature = "JAK1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["STAT1",],
                                  Feature = "STAT1")) %>%
    dplyr::rows_append(data.frame(CB = colnames(seuratObject),
                                  Embeddings(seuratObject[["umap"]]),
                                  seuratObject@meta.data,
                                  expr = GetAssayData(seuratObject, layers = "data")["IL6",],
                                  Feature = "IL6")) %>%
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
              legend.position = "bottom",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              strip.text.x = element_text(face = "bold"))
    }
)
figR1P6B <- ggarrange(plotlist = figR1P6B_list, nrow = 3, ncol = 3)

ggsave(plot = figR1P6A, filename = "figR1P6A.pdf", width = 6, height = 7, units = "in")
ggsave(plot = figR1P6B, filename = "figR1P6B.pdf", width = 9, height = 10, units = "in")
ggsave(plot = figR1P6B, filename = "figR1P6B.pdf", width = 9, height = 10.5, units = "in")
