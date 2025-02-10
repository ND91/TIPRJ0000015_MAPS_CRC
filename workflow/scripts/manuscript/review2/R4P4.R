#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 4, point 4.

BiocManager::install(c("harmony"))

require(Seurat)
require(harmony)
require(dplyr)
require(ggplot2)
require(ggrastr)

hc_pbmc_pf_seuratobject_rds <- "output/q1_pf_characterization/subsets/hc_pbmc_pf_SeuratObject.Rds"
hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seuratobject_rds <- "output/q3_pm_tx_characterization/subsets/hc_crcpmp_pbmc_pf_tx_paired_monocytes_macrophages_SeuratObject.Rds"
celltype_markers_xlsx <- "config/order/celltype_markers.xlsx"

hc_pbmc_pf_seurat <- readRDS(hc_pbmc_pf_seuratobject_rds)
hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat <- readRDS(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seuratobject_rds)
DefaultAssay(hc_pbmc_pf_seurat) <- "RNA"
DefaultAssay(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat) <- "RNA"

options(future.globals.maxSize= 20000*1024^2)

tissue_colors <- c(PF = "#224FBD", PBMC = "#F30000", TX = "#808080")
tissue_group_colors <- c("PBMC HC" = "#58508d", "PBMC CRC+" = "#bc5090", "PF HC" = "#ff6361", "PF CRC+" = "#ffa600", "TX CRC+" = "#003f5c")

celltype_order_l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(hc_pbmc_pf_seurat@meta.data[,"manual_l2"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_monocytes_macrophages <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat@meta.data[,"manual_l2"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_colors <- celltype_order_l2$color
names(manual_l2_colors) <- celltype_order_l2$celltype

manual_l2_monocytes_macrophages_colors <- manual_l2_monocytes_macrophages$color
names(manual_l2_monocytes_macrophages_colors) <- manual_l2_monocytes_macrophages$celltype

pbmc_seurat <- hc_pbmc_pf_seurat[,hc_pbmc_pf_seurat@meta.data$Tissue == "PBMC"]
pbmc_seurat <- DietSeurat(pbmc_seurat, counts = T, data = T, scale.data = F)
pbmc_seurat <- SCTransform(pbmc_seurat, conserve.memory = T)
pbmc_seurat <- RunPCA(object = pbmc_seurat, npcs = 100, seed.use = 4132)

pf_seurat <- hc_pbmc_pf_seurat[,hc_pbmc_pf_seurat@meta.data$Tissue == "PF"]
pf_seurat <- DietSeurat(pf_seurat, counts = T, data = T, scale.data = F)
pf_seurat <- SCTransform(pf_seurat, conserve.memory = T)
pf_seurat <- RunPCA(object = pf_seurat, npcs = 100, seed.use = 21352)

# PBMC and PF

## Harmony

hc_pbmc_pf_seurat_harmony <- RunHarmony(hc_pbmc_pf_seurat, "Tissue")
hc_pbmc_pf_seurat_harmony <- FindNeighbors(hc_pbmc_pf_seurat_harmony, reduction = "harmony", dims = 1:45)
hc_pbmc_pf_seurat_harmony <- RunUMAP(hc_pbmc_pf_seurat_harmony, dims = 1:45, seed.use = 398772, reduction = "harmony")

DimPlot(hc_pbmc_pf_seurat_harmony, group.by = "Tissue", reduction = "umap")
DimPlot(hc_pbmc_pf_seurat_harmony, group.by = "manual_l2", reduction = "umap", split.by = "Tissue")
DimPlot(hc_pbmc_pf_seurat_harmony, group.by = "manual_l2")

harmony_umap_hc_pbmc_pf_coltissue_ggplotobj <- data.frame(CB = colnames(hc_pbmc_pf_seurat_harmony),
                                                          Embeddings(hc_pbmc_pf_seurat_harmony[["umap"]]),
                                                          hc_pbmc_pf_seurat_harmony@meta.data) %>%
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = Tissue)) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = tissue_colors) +
  # facet_wrap(~Tissue) +
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

pdf(file.path(figR4P4Dir, "harmony_umap_hc_pbmc_pf_coltissue.pdf"), width = 7.5, height = 7.5)
print(harmony_umap_hc_pbmc_pf_coltissue_ggplotobj)
dev.off()

harmony_umap_hc_pbmc_pf_coll2_splittissue_ggplotobj <- data.frame(CB = colnames(hc_pbmc_pf_seurat_harmony),
                                                                  Embeddings(hc_pbmc_pf_seurat_harmony[["umap"]]),
                                                                  hc_pbmc_pf_seurat_harmony@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l2, levels = celltype_order_l2$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number)) %>% 
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = Tissue)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = . %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = umap_1),
                               y = median(x = umap_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   alpha = 1,
                   # size = 8,
                   show.legend = F,
                   label.size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l2_colors) +
  # scale_color_manual(values = tissue_colors) +
  facet_wrap(~Tissue) +
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

pdf(file.path(figR4P4Dir, "harmony_umap_hc_pbmc_pf_coll2_splittissue.pdf"), width = 15, height = 7.5)
print(harmony_umap_hc_pbmc_pf_coll2_splittissue_ggplotobj)
dev.off()

harmony_umap_hc_pbmc_pf_coll2_splittissue_monomacs_ggplotobj <- data.frame(CB = colnames(hc_pbmc_pf_seurat_harmony),
                                                                           Embeddings(hc_pbmc_pf_seurat_harmony[["umap"]]),
                                                                           hc_pbmc_pf_seurat_harmony@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l2, levels = celltype_order_l2$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number)) %>% 
  dplyr::filter(manual_l2 %in% c("Monocytes", "Macrophages")) %>%
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = . %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = umap_1),
                               y = median(x = umap_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   alpha = 1,
                   # size = 8,
                   show.legend = F,
                   label.size = 0.25) +
  xlim(-5, 2) +
  ylim(-12, -4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l2_colors) +
  # scale_color_manual(values = tissue_colors) +
  facet_grid(.~Tissue) +
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

pdf(file.path(figR4P4Dir, "harmony_umap_hc_pbmc_pf_coll2_splittissue_monomacs.pdf"), width = 6, height = 3.5)
print(harmony_umap_hc_pbmc_pf_coll2_splittissue_monomacs_ggplotobj)
dev.off()

umap_hc_pf_monomacrophages_colvsig4_ggplotobj <- data.frame(CB = colnames(hc_pbmc_pf_seurat_harmony),
                                                            Embeddings(hc_pbmc_pf_seurat_harmony[["umap"]]),
                                                            hc_pbmc_pf_seurat_harmony@meta.data,
                                                            expr = GetAssayData(hc_pbmc_pf_seurat_harmony, assay = "RNA")["VSIG4",],
                                                            Feature = "VSIG4",
                                                            Modality = "RNA") %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l2, levels = celltype_order_l2$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number)) %>% 
  dplyr::filter(manual_l2 %in% c("Monocytes", "Macrophages")) %>%
  dplyr::arrange(expr) %>%
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = . %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T) +
  geom_label_repel(data = . %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = umap_1),
                               y = median(x = umap_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   alpha = 1,
                   # size = 8,
                   show.legend = F,
                   label.size = 0.25) +
  xlim(-5, 2) +
  ylim(-12, -4) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_color_viridis() +
  facet_grid(.~Tissue) +
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

pdf(file.path(figR4P4Dir, "harmony_umap_hc_pf_monomacrophages_colvsig4.pdf"), width = 6, height = 3)
print(umap_hc_pf_monomacrophages_colvsig4_ggplotobj)
dev.off()

# PBMC PF and TX

hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony <- RunHarmony(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat, "Tissue")
hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony <- FindNeighbors(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony, reduction = "harmony", dims = 1:46)
hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony <- RunUMAP(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony, dims = 1:46, seed.use = 7869243, reduction = "harmony")

umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coll2_ggplotobj <- data.frame(
  CB = colnames(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony),
  Embeddings(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony[["umap"]]),
  hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony@meta.data) %>%
  ggplot(aes(x = umap_1, y = umap_2, col = manual_l2)) +
  geom_point_rast(show.legend = T, size = 1) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l2_monocytes_macrophages_colors) +
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

pdf(file.path(figR4P4Dir, "umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coll2.pdf"), width = 7.5, height = 7.5)
print(umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coll2_ggplotobj)
dev.off()

umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coltissuegroup_ggplotobj <- data.frame(
  CB = colnames(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony),
  Embeddings(hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony[["umap"]]),
  hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_seurat_harmony@meta.data) %>%
  dplyr::mutate(Tissue_group = factor(paste0(Tissue, " ", Group), levels = c("PBMC HC", "PBMC CRC+", "PF HC", "PF CRC+", "TX CRC+"))) %>%
  ggplot(aes(x = umap_1, y = umap_2, col = Tissue_group)) +
  geom_point_rast(show.legend = T, size = 1) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = tissue_group_colors) +
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

pdf(file.path(figR4P4Dir, "umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coltissuegroup.pdf"), width = 7.5, height = 7.5)
print(umap_hc_crcpmp_pbmc_pf_tx_monocytes_macrophages_coltissuegroup_ggplotobj)
dev.off()