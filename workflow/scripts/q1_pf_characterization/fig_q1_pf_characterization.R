#!/usr/bin/env Rscript

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(ggpubr)
require(viridis)

seurat_pf_rds <- args[1]
seurat_pf_t_rds <- args[2]
seurat_pf_nkilc_rds <- args[3]
seurat_pf_myeloid_rds <- args[4]
seurat_pbmc_pf_rds <- args[5]
seurat_pbmc_pf_liver_colon_rds <- args[6]
pfvpbmc_degs_l3_rds <- args[7]
pfvpbmc_fgsea_l3_rds <- args[8]
livervpf_dacs_l3_csv <- args[9]
colonvpf_dacs_l3_csv <- args[10]
celltype_markers_xlsx <- args[11]
cite_markers_xlsx <- args[12]
dacs_pbmcvpf_l3rl0_rds <- args[13]

seurat_pf <- readRDS(seurat_pf_rds)
seurat_pf_t <- readRDS(seurat_pf_t_rds)
seurat_pf_nkilc <- readRDS(seurat_pf_nkilc_rds)
seurat_pf_myeloid <- readRDS(seurat_pf_myeloid_rds)
seurat_pbmc_pf <- readRDS(seurat_pbmc_pf_rds)
seurat_pbmc_pf_liver_colon <- readRDS(seurat_pbmc_pf_liver_colon_rds)
pfvpbmc_degs_l3 <- readRDS(pfvpbmc_degs_l3_rds)
pfvpbmc_fgsea_l3 <- readRDS(pfvpbmc_fgsea_l3_rds)
livervpf_dacs_l3 <- read.csv(livervpf_dacs_l3_csv)
colonvpf_dacs_l3 <- read.csv(colonvpf_dacs_l3_csv)

dacs_pbmcvpf_l3rl0 <- readRDS(dacs_pbmcvpf_l3rl0_rds)

manual_l1_order <- readxl::read_excel(celltype_markers_xlsx) %>%
  dplyr::filter(level == "manual_l1") %>%
  dplyr::pull(celltype) %>%
  unique()

manual_l3_order <- readxl::read_excel(celltype_markers_xlsx) %>%
  dplyr::filter(level == "manual_l3") %>%
  dplyr::pull(celltype) %>%
  unique()

tissue_colors <- c(PF = "#224FBD", PBMC = "#F30000", Colon = "#7B3F00", Liver = "#AA336A")

celltype_order_pf_l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_pf@meta.data[,"manual_l2"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_colors <- celltype_order_pf_l2$color
names(manual_l2_colors) <- celltype_order_pf_l2$celltype

celltype_order_pf_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_pf@meta.data[,"manual_l3"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))


celltype_order_l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(!is.na(color)) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_colors_full <- celltype_order_l2$color
names(manual_l2_colors_full) <- celltype_order_l2$celltype

# Fig UMAP HC PF color: manual_l1

celltype_order_pf_l1 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l1",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_pf@meta.data[,"manual_l1"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

celltype_colors_pf_l1_list <- celltype_order_pf_l1$color
names(celltype_colors_pf_l1_list) <- celltype_order_pf_l1$celltype

umap_hc_pf_l1_df <- data.frame(CB = colnames(seurat_pf),
                               Embeddings(seurat_pf[["wnn.umap"]]),
                               seurat_pf@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l1, levels = celltype_order_pf_l1$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l1$number_subset))

umap_hc_pf_l1_ggplotobj <- umap_hc_pf_l1_df %>% 
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "black") +
  geom_point_rast(show.legend = F, size = 0.25, aes(col = celltype)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_pf_l1_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = wnnUMAP_1),
                               y = median(x = wnnUMAP_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_colors_pf_l1_list) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_l1_together.pdf", width = 7.5, height = 7.5)
print(umap_hc_pf_l1_ggplotobj)
dev.off()

# Fig UMAP HC PF color: manual_l3

umap_hc_pf_together_l3_df <- data.frame(CB = colnames(seurat_pf),
                                        Embeddings(seurat_pf[["wnn.umap"]]),
                                        seurat_pf@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

number_l3_pf_list <- celltype_order_pf_l3$color
names(number_l3_pf_list) <- celltype_order_pf_l3$number_subset

umap_hc_pf_l3_together_ggplotobj <- umap_hc_pf_together_l3_df %>% 
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_text_repel(data = umap_hc_pf_together_l3_df %>%
                    dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                    summarize(x = median(x = wnnUMAP_1),
                              y = median(x = wnnUMAP_2)),
                  mapping = aes(label = celltype_number, x = x, y = y),
                  alpha = 1,
                  size = 6,
                  show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = number_l3_pf_list) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_l3_together.pdf", width = 11, height = 10)
print(umap_hc_pf_l3_together_ggplotobj)
dev.off()

# Fig UMAP HC PF color: manual_l3, split: manual_l1

celltype_order_pf_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_pf@meta.data[,"manual_l3"])) %>%
  dplyr::mutate(number = 1:nrow(.),
                number_subset = paste0(1:nrow(.), ". ", celltype))

celltype_number_colors_pf_l3_list <- celltype_order_pf_l3$color
names(celltype_number_colors_pf_l3_list) <- celltype_order_pf_l3$number_subset

celltype_colors_pf_l3_list <- celltype_order_pf_l3$color
names(celltype_colors_pf_l3_list) <- celltype_order_pf_l3$celltype

umap_hc_pf_l3_df <- data.frame(CB = colnames(seurat_pf_t),
                               Embeddings(seurat_pf_t[["wnn.umap"]]),
                               seurat_pf_t@meta.data) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_nkilc),
                                Embeddings(seurat_pf_nkilc[["wnn.umap"]]),
                                seurat_pf_nkilc@meta.data)) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data)) %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_t_gplotobj <- umap_hc_pf_l3_df %>% 
  dplyr::filter(manual_l1 == "T") %>%
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_pf_l3_df %>%
                     dplyr::filter(manual_l1 == "T") %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = wnnUMAP_1),
                               y = median(x = wnnUMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_number_colors_pf_l3_list) +
  facet_wrap(~manual_l1, nrow = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

umap_hc_pf_nkilc_gplotobj <- umap_hc_pf_l3_df %>% 
  dplyr::filter(manual_l1 == "NK/ILC") %>%
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_pf_l3_df %>%
                     dplyr::filter(manual_l1 == "NK/ILC") %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = wnnUMAP_1),
                               y = median(x = wnnUMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_number_colors_pf_l3_list) +
  facet_wrap(~manual_l1, nrow = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

umap_hc_pf_myeloid_ggplotobj <- umap_hc_pf_l3_df %>% 
  dplyr::filter(manual_l1 == "Myeloid") %>%
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_pf_l3_df %>%
                     dplyr::filter(manual_l1 == "Myeloid") %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = wnnUMAP_1),
                               y = median(x = wnnUMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_number_colors_pf_l3_list) +
  facet_wrap(~manual_l1, nrow = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

umap_hc_pf_l3_gplotobj <- ggarrange(umap_hc_pf_t_gplotobj, 
                                    umap_hc_pf_nkilc_gplotobj,
                                    umap_hc_pf_myeloid_gplotobj,
                                    nrow = 3, 
                                    align = "hv")

# Fig UMAP HC PF color by CD3 (T), CD56 (NK/ILC), CD20 (B), and CD11b (Myeloid) at both GEX and PEX 

umap_hc_pf_cd3_cd56_cd20_cd11b_df <- data.frame(CB = colnames(seurat_pf),
                                                Embeddings(seurat_pf[["wnn.umap"]]),
                                                seurat_pf@meta.data,
                                                expr = GetAssayData(seurat_pf, assay = "RNA")["CD3D",],
                                                Feature = "CD3/CD3D",
                                                Modality = "RNA") %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "CITE")["Hu.CD3-UCHT1",],
                                Feature = "CD3/CD3D",
                                Modality = "CITE"))  %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "RNA")["NCAM1",],
                                Feature = "CD56/NCAM1",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "CITE")["Hu.CD56",],
                                Feature = "CD56/NCAM1",
                                Modality = "CITE")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "RNA")["MS4A1",],
                                Feature = "CD20/MS4A1",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "CITE")["Hu.CD20-2H7",],
                                Feature = "CD20/MS4A1",
                                Modality = "CITE")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "RNA")["ITGAM",],
                                Feature = "CD11b/ITGAM",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf),
                                Embeddings(seurat_pf[["wnn.umap"]]),
                                seurat_pf@meta.data,
                                expr = GetAssayData(seurat_pf, assay = "CITE")["Hu.CD11b",],
                                Feature = "CD11b/ITGAM",
                                Modality = "CITE")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                Feature = factor(Feature, levels = c("CD3/CD3D", "CD56/NCAM1", "CD20/MS4A1", "CD11b/ITGAM")),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_cd3_cd56_cd20_cd11b_ggplotobj <- ggplot(umap_hc_pf_cd3_cd56_cd20_cd11b_df, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = umap_hc_pf_cd3_cd56_cd20_cd11b_df %>%
                    dplyr::filter(expr>0), aes(col = expr, order = expr_rank), show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_grid(Modality~Feature, switch = "y") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_cd3_cd56_cd20_cd11b.pdf", width = 10, height=6)
print(umap_hc_pf_cd3_cd56_cd20_cd11b_ggplotobj)
dev.off()

# Fig heatmap PEX HC PF manual_l3

marker_proteins <- readxl::read_excel(cite_markers_xlsx) %>%
  dplyr::filter(Assay == "CITE")

seurat_pf@meta.data$manual_l3 <- factor(seurat_pf@meta.data$manual_l3, levels = manual_l3_order)

Idents(seurat_pf) <- "manual_l3"
seurat_pf_avexpr <- AverageExpression(seurat_pf, return.seurat = T)

marker_avepexpr_pex <- GetAssayData(seurat_pf_avexpr, assay = "CITE")[marker_proteins$FeatureID, ]

heatmap_pex_complexheatmapobj <- ComplexHeatmap::pheatmap(marker_avepexpr_pex, 
                                                          color = PurpleAndYellow(1000),
                                                          labels_row = marker_proteins$protein,
                                                          cluster_rows = F, 
                                                          cluster_cols = F, 
                                                          scale = "row",
                                                          name = "Scaled median")

pdf("heatmap_pex.pdf", width = 6.5, height=6.25)
print(heatmap_pex_complexheatmapobj)
dev.off()

# Fig dotplot GEX HC PF manual_l3 (all)

marker_genes <- readxl::read_excel(cite_markers_xlsx) %>%
  dplyr::filter(Assay == "RNA")
marker_gex <- GetAssayData(seurat_pf, assay = "RNA")[which(rownames(GetAssayData(seurat_pf, assay = "RNA")) %in% unique(marker_genes$FeatureID)), ]

marker_median_gex <- data.frame(FeatureID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-c(FeatureID), names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seurat_pf@meta.data), 
                               Celltype = seurat_pf@meta.data[,"manual_l3"]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                FeatureID = factor(FeatureID)) %>%
  dplyr::group_by(Celltype, FeatureID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, manual_l3_order),
                FeatureID = factor(FeatureID, rev(marker_genes$FeatureID)))

dotplot_gex_ggplotobj <- marker_median_gex %>% 
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = Celltype, y = FeatureID, col = Median)) +
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

pdf("dotplot_gex_full.pdf", width = 6, height=15)
print(dotplot_gex_ggplotobj)
dev.off()

# Fig dotplot GEX HC PF manual_l3 (overlapping with pex)

dotplot_gex_overlapping_ggplotobj <- marker_median_gex %>% 
  dplyr::filter(FeatureID %in% marker_proteins$gene) %>%
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = Celltype, y = FeatureID, col = Median)) +
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

pdf("dotplot_gex.pdf", width = 5.75, height=7.5)
print(dotplot_gex_overlapping_ggplotobj)
dev.off()

# Fig UMAP HC PF Myeloid color by CD163, MARCO, VSIG4, C5AR1 

umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1_df <- data.frame(CB = colnames(seurat_pf_myeloid),
                                                            Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                                            seurat_pf_myeloid@meta.data,
                                                            expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["CD163",],
                                                            Feature = "CD163",
                                                            Modality = "RNA") %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["MARCO",],
                                Feature = "MARCO",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["C5AR1",],
                                Feature = "C5AR1",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["VSIG4",],
                                Feature = "VSIG4",
                                Modality = "RNA")) %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["CD1C",],
                                Feature = "CD1C",
                                Modality = "RNA")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                Feature = factor(Feature, levels = c("CD163", "MARCO", "C5AR1", "VSIG4", "CD1C")),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1_ggplotobj <- ggplot(umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1_df, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1_df %>%
                    dplyr::filter(expr>0), aes(col = expr, order = expr_rank), show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~Feature, nrow = 1, ncol = 6) +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1.pdf", width = 18, height=3)
print(umap_hc_pf_myeloid_cd163_marco_vsig4_c5ar1_ggplotobj)
dev.off()

# Fig UMAP HC PF Myeloid color by CD163 at both GEX and PEX

umap_hc_pf_myeloid_cd163_df <- data.frame(CB = colnames(seurat_pf_myeloid),
                                          Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                          seurat_pf_myeloid@meta.data,
                                          expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["CD163",],
                                          Modality = "RNA") %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "CITE")["Hu.CD163",],
                                Modality = "CITE")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_myeloid_cd163_ggplotobj <- ggplot(umap_hc_pf_myeloid_cd163_df, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = umap_hc_pf_myeloid_cd163_df %>%
                    dplyr::filter(expr>0), aes(col = expr, order = expr_rank), show.legend = T, size = 0.25) +
  # geom_label_repel(data = umap_hc_pf_myeloid_cd163_df %>%
  #                    dplyr::group_by(manual_l2) %>%
  #                    summarize(x = median(x = wnnUMAP_1),
  #                              y = median(x = wnnUMAP_2)),
  #                  mapping = aes(label = manual_l2, x = x, y = y),
  #                  alpha = 1, 
  #                  show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~Modality, ncol = 2) +
  labs(title = "CD163 expression") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_myeloid_cd163.pdf", width = 6, height=3.5)
print(umap_hc_pf_myeloid_cd163_ggplotobj)
dev.off()

# Fig UMAP HC PF Myeloid color by C1QA and VCAN

umap_hc_pf_myeloid_c1qa_vcan_df <- data.frame(CB = colnames(seurat_pf_myeloid),
                                              Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                              seurat_pf_myeloid@meta.data,
                                              expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["C1QA",],
                                              Gene = "C1QA") %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_pf_myeloid),
                                Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                seurat_pf_myeloid@meta.data,
                                expr = GetAssayData(seurat_pf_myeloid, assay = "RNA")["VCAN",],
                                Gene = "VCAN")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_myeloid_c1qa_vcan_ggplotobj <- ggplot(umap_hc_pf_myeloid_c1qa_vcan_df, aes(x = wnnUMAP_1, y = wnnUMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = umap_hc_pf_myeloid_c1qa_vcan_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.25) +
  facet_wrap(~Gene, ncol = 2) +
  #labs(title = "Gene expression") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf("umap_hc_pf_myeloid_c1qa_vcan.pdf", width = 6, height=3.5)
print(umap_hc_pf_myeloid_c1qa_vcan_ggplotobj)
dev.off()

# Fig UMAP HC PF Macrophages color by C1QA and VCAN

densityplot_hc_pf_macrophages_c1qa_vcan_df <- data.frame(CB = colnames(seurat_pf_myeloid),
                                                         Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                                         seurat_pf_myeloid@meta.data,
                                                         VCAN = GetAssayData(seurat_pf_myeloid, assay = "RNA")["VCAN",],
                                                         C1QA = GetAssayData(seurat_pf_myeloid, assay = "RNA")["C1QA",]) %>%
  dplyr::filter(manual_l2 == "Macrophages")

densityplot_hc_pf_macrophages_c1qa_vcan_ggplotobj <- densityplot_hc_pf_macrophages_c1qa_vcan_df %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset)) %>%
  ggplot(aes(x = VCAN, y = C1QA)) +
  #geom_point_rast(alpha = 0.15) +
  geom_density_2d_filled(alpha = 0.9) +
  labs(title = "Macrophages: Overlap C1QA and VCAN expression") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

pdf("densityplot_hc_pf_macrophages_c1qa_vcan.pdf", width = 5, height=5)
print(densityplot_hc_pf_macrophages_c1qa_vcan_ggplotobj)
dev.off()

# Fig UMAP HC PF color by manual_l3

umap_hc_pf_myeloid_l3_df <- data.frame(CB = colnames(seurat_pf_myeloid),
                                       Embeddings(seurat_pf_myeloid[["wnn.umap"]]),
                                       seurat_pf_myeloid@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_pf_myeloid_l3_ggplotobj <- umap_hc_pf_myeloid_l3_df %>% 
  ggplot(aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_pf_myeloid_l3_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = wnnUMAP_1),
                               y = median(x = wnnUMAP_2)),
                   mapping = aes(label = celltype, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = celltype_colors_pf_l3_list) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

# Fig boxplot manual_l2 proportions relative to CD45+

l2rl0_proportions_df <- seurat_pf@meta.data %>%
  dplyr::group_by(manual_l2, SampleID, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(SampleID), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) 

l2rl0_proportions_ggplotobj <- l2rl0_proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l2, -Ncellperc), y = Ncellperc, fill = manual_l2)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point_rast(aes(shape = Donor), size = 3) +
  scale_fill_manual(values = manual_l2_colors) +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("l2rl0_proportions.pdf", width = 7.5, height=5)
print(l2rl0_proportions_ggplotobj)
dev.off()

# Fig UMAP HC PBMC PF color by tissue

umap_hc_pbmc_pf <- data.frame(CB = colnames(seurat_pbmc_pf),
                              Embeddings(seurat_pbmc_pf[["umap"]]),
                              seurat_pbmc_pf@meta.data)

umap_coltissue_ggplotobj <- ggplot(umap_hc_pbmc_pf, aes(x = UMAP_1, y = UMAP_2, col = Tissue, partition = Tissue)) +
  #geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  #geom_point_rast(show.legend = T, size = 0.25) * (blend("lighten") + blend("multiply", alpha = 0.5))+
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = tissue_colors) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

cairo_pdf("umap_coltissue.pdf", width = 7.5, height=7.5)
print(umap_coltissue_ggplotobj)
dev.off()

umap_coltissue_splittissue_ggplotobj <- ggplot(umap_hc_pbmc_pf, aes(x = UMAP_1, y = UMAP_2, col = Tissue, partition = Tissue)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = tissue_colors) +
  facet_wrap(~Tissue, nrow = 1) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

cairo_pdf("umap_coltissue_splittissue.pdf", width = 15, height=7.5)
print(umap_coltissue_splittissue_ggplotobj)
dev.off()

# Fig boxplot PBMC PF manual_l3 relative l0

pbmc_pf_l3rl0_proportions_df <- seurat_pbmc_pf@meta.data %>%
  dplyr::group_by(manual_l3, SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(unique(seurat_pbmc_pf@meta.data[,c("manual_l1", "manual_l3")]), by = "manual_l3") %>%
  dplyr::mutate(manual_l1 = factor(manual_l1, levels = manual_l1_order)) %>%
  dplyr::left_join(dacs_pbmcvpf_l3rl0[[1]]$da, by = c("manual_l3" = "BaselineProp.clusters"))

boxplot_pbmc_pf_l3rl0_together_ggplotobj <- pbmc_pf_l3rl0_proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellperc), y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_grid(.~ manual_l1, scales = "free", space="free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

boxplot_pbmc_pf_l3rl0_together_ggplotobj <- ggarrange(plotlist = lapply(split(pbmc_pf_l3rl0_proportions_df, factor(pbmc_pf_l3rl0_proportions_df$manual_l1, levels = c("T", "NK/ILC", "B", "Myeloid"))), function(manual_l1){
  manual_l1 %>%
    ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellperc), y = Ncellperc, col = Tissue)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
    facet_wrap(.~ manual_l1) +
    labs(y = "%CD45+ cells") +
    theme_bw() +
    scale_color_manual(values = tissue_colors) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.pos = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}), widths = c(15, 4, 4, 8), nrow = 1, align = "hv", legend = "bottom", common.legend = T)

pdf("boxplot_pbmc_pf_l3rl0_together.pdf", width = 12.5, height = 5)
print(boxplot_pbmc_pf_l3rl0_together_ggplotobj)
dev.off()

boxplot_pbmc_pf_l3rl0_patanno_separate_ggplotobj <- pbmc_pf_l3rl0_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l3, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label))) %>%
  ggplot(aes(x = Tissue, y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_pbmc_pf_l3rl0_separate.pdf", width = 10, height = 10)
print(boxplot_pbmc_pf_l3rl0_patanno_separate_ggplotobj)
dev.off()

boxplot_pbmc_pf_patanno_together_ggplotobj <- pbmc_pf_l3rl0_proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellperc), y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
  facet_grid(.~ manual_l1, scales = "free", space='free') +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_pbmc_pf_patanno_together.pdf", width = 12.5, height = 5)
print(boxplot_pbmc_pf_patanno_together_ggplotobj)
dev.off()

boxplot_pbmc_pf_patanno_separate_ggplotobj <- l3rl0_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l3, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label))) %>%
  ggplot(aes(x = Tissue, y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5, aes(shape = Donor)) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_pbmc_pf_patanno_separate.pdf", width = 10, height = 10)
print(boxplot_pbmc_pf_patanno_separate_ggplotobj)
dev.off()

dpws <- sort(table(unlist(lapply(pfvpbmc_fgsea_l3, function(celltype){
  celltype %>%
    dplyr::select(1:7) %>%
    data.frame() %>%
    dplyr::filter(padj<0.05) %>% 
    dplyr::pull(pathway)
}))))

# Fig stackedbarplot PBMC, PF, Liver and Colon

proportions_l2rl0_tissue_df <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::group_by(manual_l2, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Tissue) %>%
  dplyr::mutate(manual_l2 = factor(manual_l2, levels = celltype_order_l2$celltype),
                Nsample = sum(Ncells),
                Ncellprop = Ncells/Nsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver"))) %>%
  dplyr::rename(Celltype = manual_l2)

stackedbarplot_l2rl0_tissue_ggplotobj <- proportions_l2rl0_tissue_df %>%
  ggplot(aes(x = Tissue, y = Ncellperc)) +
  geom_bar(position="stack", stat="identity", aes(fill = Celltype)) +
  labs(y = "%CD45+") +
  scale_fill_manual(values = manual_l2_colors_full) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("stackedbarplot_l2rl0_tissue.pdf", width = 7.5, height = 7.5)
print(stackedbarplot_l2rl0_tissue_ggplotobj)
dev.off()

# Fig boxplot l2 PBMC, PF, Liver and Colon 

proportions_l2rl0_tissue_df <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::group_by(manual_l2, Donor, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, Donor), tidyr::nesting(manual_l2), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Donor) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", Donor),
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver")))

boxplot_l2rl0_tissue_separate_ggplotobj <- proportions_l2rl0_tissue_df %>%
  ggplot(aes(x = Tissue, y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~manual_l2, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_l2rl0_tissue_separate.pdf", width = 10, height = 10)
print(boxplot_l2rl0_tissue_separate_ggplotobj)
dev.off()

# Fig scatterplot l3 liverpf and colonvpf

pfvlivervcolon_df <- livervpf_dacs_l3 %>%
  dplyr::mutate(Tstatistic=Tstatistic*-1) %>%
  dplyr::select(BaselineProp.clusters, Tstatistic, P.Value, FDR) %>%
  dplyr::rename(Celltype = BaselineProp.clusters,
                tstat_livervpf = Tstatistic, 
                pvalue_livervpf = P.Value,
                fdr_livervpf = FDR) %>%
  dplyr::inner_join(colonvpf_dacs_l3 %>%
                     dplyr::mutate(Tstatistic=Tstatistic*-1) %>%
                     dplyr::select(BaselineProp.clusters, Tstatistic, P.Value, FDR) %>%
                     dplyr::rename(Celltype = BaselineProp.clusters,
                                   tstat_colonvpf = Tstatistic, 
                                   pvalue_colonvpf = P.Value,
                                   fdr_colonvpf = FDR), by = "Celltype")
scatterplot_pfvlivervcolon_ggplotobj <- pfvlivervcolon_df %>%
  dplyr::mutate(Significant = ifelse(fdr_livervpf<0.05 & fdr_colonvpf<0.05, "Significant", "NS"),
                label = ifelse(Significant == "Significant", Celltype, NA)) %>%
  ggplot(aes(x = tstat_colonvpf, y = tstat_livervpf)) +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_point(aes(alpha = Significant)) +
  geom_label_repel(aes(label = label)) +
  labs(x = "Colon vs PF",
       y = "Liver vs PF") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

pdf("scatterplot_pfvlivervcolon.pdf", width = 7.5, height = 7.5)
print(scatterplot_pfvlivervcolon_ggplotobj)
dev.off()

# Fig boxplot l3 PBMC, PF, Liver and Colon 

proportions_l3rl0_tissue_df <- seurat_pbmc_pf_liver_colon@meta.data %>%
  dplyr::filter(Tissue %in% c("PF", "Colon", "Liver")) %>%
  dplyr::group_by(manual_l3, Donor, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, Donor), tidyr::nesting(manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(Donor) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", Donor),
                Tissue = factor(Tissue, levels = c("PBMC", "PF", "Colon", "Liver")))

boxplot_l3rl0_tissue_separate_ggplotobj <- proportions_l3rl0_tissue_df %>%
  ggplot(aes(x = Tissue, y = Ncellperc, col = Tissue)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~manual_l3, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = tissue_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_l3rl0_tissue_separate.pdf", width = 15, height = 15)
print(boxplot_l3rl0_tissue_separate_ggplotobj)
dev.off()

hc_pf_l3rl1_proportions_df <- seurat_pf@meta.data %>%
  dplyr::group_by(manual_l3, manual_l1, SampleID, Tissue, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Tissue, SampleID), tidyr::nesting(manual_l3, manual_l1), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID, manual_l1) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID))
