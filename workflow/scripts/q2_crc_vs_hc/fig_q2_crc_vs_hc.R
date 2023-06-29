#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of of the abundances relative to the reference level.

require(Seurat)
require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)
require(ggpubr)
require(viridis)
require(ggblend)

seurat_hc_crcpmp_pf_rds <- args[1]
seurat_hc_crcpmp_pf_macrophages_rds <- args[2]
dacs_pf_l3rl0_crcpmpvhc_rds <- args[3]
dacs_pf_l4rl0_crcpmpvhc_rds <- args[4]
degs_pf_l3_crcpmpvhc_rds <- args[5]
degs_pf_l4_crcpmpvhc_rds <- args[6]
fgsea_pf_l4_crcpmpvhc_rds <- args[7]
celltype_markers_xlsx <- args[8]

macrophage_l4 <- c("Macrophages VCAN+", "Macrophages C1Q+", "Macrophages C1Q+SPP1+", "Macrophages C1Q+MAFbright", "Macrophages VCAN+C1Q+")

seurat_hc_crcpmp_pf <- readRDS(seurat_hc_crcpmp_pf_rds)
seurat_hc_crcpmp_pf_macrophages <- readRDS(seurat_hc_crcpmp_pf_macrophages_rds)
dacs_pf_l3rl0_crcpmpvhc <- readRDS(dacs_pf_l3rl0_crcpmpvhc_rds)
dacs_pf_l4rl0_crcpmpvhc <- readRDS(dacs_pf_l4rl0_crcpmpvhc_rds)
degs_pf_l3_crcpmpvhc <- readRDS(degs_pf_l3_crcpmpvhc_rds)
degs_pf_l4_crcpmpvhc <- readRDS(degs_pf_l4_crcpmpvhc_rds)
fgsea_pf_l4_crcpmpvhc <- readRDS(fgsea_pf_l4_crcpmpvhc_rds)

tissue_colors <- c(PF = "#224FBD", PBMC = "#F30000")
group_colors <- c(`HC` = "#93CEC1", `CRC+` = "#996633")

manual_l1_order <- readxl::read_excel(celltype_markers_xlsx) %>%
  dplyr::filter(level == "manual_l1") %>%
  dplyr::pull(celltype) %>%
  unique()

manual_l3_order <- readxl::read_excel(celltype_markers_xlsx) %>%
  dplyr::filter(level == "manual_l3") %>%
  dplyr::pull(celltype) %>%
  unique()


celltype_order_pf_l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_hc_crcpmp_pf@meta.data[,"manual_l2"])) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_colors <- celltype_order_pf_l2$color
names(manual_l2_colors) <- celltype_order_pf_l2$celltype

celltype_order_l2 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l2",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(!is.na(color)) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

manual_l2_colors_full <- celltype_order_l2$color
names(manual_l2_colors_full) <- celltype_order_l2$celltype

celltype_order_pf_l3 <- readxl::read_excel(celltype_markers_xlsx, col_names = T) %>%
  dplyr::filter(level == "manual_l3",
                modality == "gene") %>%
  dplyr::select(celltype, color) %>%
  unique() %>%
  dplyr::filter(celltype %in% unique(seurat_hc_crcpmp_pf@meta.data[,"manual_l3"]),
                !is.na(color)) %>%
  dplyr::mutate(number_subset = paste0(1:nrow(.), ". ", celltype))

celltype_colors_pf_l3_list <- celltype_order_pf_l3$color
names(celltype_colors_pf_l3_list) <- celltype_order_pf_l3$number_subset

# Fig UMAP HC CRC PM+ PF color by manual_l3

umap_hc_crcpmp_pf_l3_df <- data.frame(CB = colnames(seurat_hc_crcpmp_pf),
                                      Embeddings(seurat_hc_crcpmp_pf[["umap"]]),
                                      seurat_hc_crcpmp_pf@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_crcpmp_pf_l3_plotobj <- umap_hc_crcpmp_pf_l3_df %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_crcpmp_pf_l3_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
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

pdf("umap_hc_crcpmp_pf_l3.pdf", width = 8, height=9)
print(umap_hc_crcpmp_pf_l3_plotobj)
dev.off()

# Fig UMAP HC CRC PM+ PF color by manual_l3

umap_hc_crcpmp_pf_macrophages_l4_df <- data.frame(CB = colnames(seurat_hc_crcpmp_pf_macrophages),
                                                  Embeddings(seurat_hc_crcpmp_pf_macrophages[["umap"]]),
                                                  seurat_hc_crcpmp_pf_macrophages@meta.data) %>%
  dplyr::mutate(celltype = factor(manual_l4, levels = macrophage_l4),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number))

umap_hc_crcpmp_pf_macrophages_l4_ggplotobj <- umap_hc_crcpmp_pf_macrophages_l4_df %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_crcpmp_pf_macrophages_l4_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  #scale_color_manual(values = celltype_colors_pf_l3_list) +
  facet_wrap(~Group, nrow = 2) +
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

pdf("umap_hc_crcpmp_pf_macrophages_l4.pdf", width = 5, height=10.5)
print(umap_hc_crcpmp_pf_macrophages_l4_ggplotobj)
dev.off()

# Fig UMAP HC CRC PM+ PF color by manual_l4

umap_hc_crcpmp_pf_macrophages_colgroup_ggplotobj <- umap_hc_crcpmp_pf_macrophages_l4_df %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, col = Group)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25, aes(col = celltype_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = umap_hc_crcpmp_pf_macrophages_l4_df %>%
                     dplyr::group_by(celltype, celltype_number, celltype_w_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = celltype_number, x = x, y = y),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  #scale_color_manual(values = celltype_colors_pf_l3_list) +
  facet_wrap(~Group, nrow = 2) +
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

pdf("umap_hc_crcpmp_pf_macrophages_l4.pdf", width = 5, height=10.5)
print(umap_hc_crcpmp_pf_macrophages_l4_ggplotobj)
dev.off()

# Fig UMAP HC CRC PM+ PF color by Donor

umap_hc_crcpmp_pf_macrophages_coldonor_ggplotobj <- ggplot(umap_hc_crcpmp_pf_macrophages_l4_df, aes(x = UMAP_1, y = UMAP_2, col = Donor, partition = Donor)) +
  geom_point_rast(show.legend = T, size = 0.25) * (blend("lighten") + blend("multiply", alpha = 0.5))+
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  #scale_color_manual(values = group_colors) +
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

cairo_pdf("umap_hc_crcpmp_pf_macrophages_coldonor.pdf", width = 7.5, height=7.5)
print(umap_hc_crcpmp_pf_macrophages_coldonor_ggplotobj)
dev.off()

# Fig Heatmap HC CRC PM+ PF macrophages manaul_l4 by donor

pdf("heatmap_hc_crcpmp_pf_macrophages_l4_proportion_donor.pdf", width = 3.5, height=10)
pheatmap::pheatmap(prop.table(table(umap_hc_crcpmp_pf_macrophages_l4_df$Donor, umap_hc_crcpmp_pf_macrophages_l4_df$manual_l4), margin = 2), display_numbers = T)
dev.off()

# Fig UMAP HC CRC PM+ PF color by C1QA, VCAN, MAF, SPP1

umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1_df <- data.frame(CB = colnames(seurat_hc_crcpmp_pf_macrophages),
                                                                  Embeddings(seurat_hc_crcpmp_pf_macrophages[["umap"]]),
                                                                  seurat_hc_crcpmp_pf_macrophages@meta.data,
                                                                  expr = GetAssayData(seurat_hc_crcpmp_pf_macrophages, assay = "RNA")["C1QA",],
                                                                  Gene = "C1QA") %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_hc_crcpmp_pf_macrophages),
                                Embeddings(seurat_hc_crcpmp_pf_macrophages[["umap"]]),
                                seurat_hc_crcpmp_pf_macrophages@meta.data,
                                expr = GetAssayData(seurat_hc_crcpmp_pf_macrophages, assay = "RNA")["VCAN",],
                                Gene = "VCAN"))  %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_hc_crcpmp_pf_macrophages),
                                Embeddings(seurat_hc_crcpmp_pf_macrophages[["umap"]]),
                                seurat_hc_crcpmp_pf_macrophages@meta.data,
                                expr = GetAssayData(seurat_hc_crcpmp_pf_macrophages, assay = "RNA")["MAF",],
                                Gene = "MAF"))  %>%
  dplyr::rows_append(data.frame(CB = colnames(seurat_hc_crcpmp_pf_macrophages),
                                Embeddings(seurat_hc_crcpmp_pf_macrophages[["umap"]]),
                                seurat_hc_crcpmp_pf_macrophages@meta.data,
                                expr = GetAssayData(seurat_hc_crcpmp_pf_macrophages, assay = "RNA")["SPP1",],
                                Gene = "SPP1")) %>%
  dplyr::mutate(expr_rank = rank(expr, ties.method="first"),
                celltype = factor(manual_l3, levels = celltype_order_pf_l3$celltype),
                celltype_number = as.numeric(celltype),
                celltype_w_number = paste0(as.numeric(celltype), ". ", celltype),
                celltype_w_number = factor(celltype_w_number, levels = celltype_order_pf_l3$number_subset))

umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1_ggplotobj <- ggplot(umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1_df %>%
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

pdf("umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1.pdf", width = 6, height=6.5)
print(umap_hc_crcpmp_pf_macrophages_c1qa_vcan_maf_spp1_ggplotobj)
dev.off()

# Fig boxplot CRC PM+ vs HC PF manual_l3 relative l0

crcpmpvhc_pf_l3rl0_proportions_df <- seurat_hc_crcpmp_pf@meta.data %>%
  dplyr::group_by(manual_l3, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(unique(seurat_hc_crcpmp_pf@meta.data[,c("manual_l1", "manual_l3")]), by = "manual_l3") %>%
  dplyr::mutate(manual_l1 = factor(manual_l1, levels = manual_l1_order)) %>%
  dplyr::left_join(dacs_pf_l3rl0_crcpmpvhc[[1]]$da, by = c("manual_l3" = "BaselineProp.clusters"))

boxplot_crcpmpvhc_pf_together_ggplotobj <- crcpmpvhc_pf_l3rl0_proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellperc), y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_grid(.~ manual_l1, scales = "free", space="free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_together.pdf", width = 12.5, height = 5)
print(boxplot_crcpmpvhc_pf_together_ggplotobj)
dev.off()

boxplot_crcpmpvhc_pf_separate_ggplotobj <- crcpmpvhc_pf_l3rl0_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l3, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label)),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  ggplot(aes(x = Group, y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_l3rl0_separate.pdf", width = 10, height = 12.5)
print(boxplot_crcpmpvhc_pf_separate_ggplotobj)
dev.off()

# Fig boxplot CRC PM+ vs HC PF manual_l4 relative l0

crcpmpvhc_pf_l4rl0_proportions_df <- seurat_hc_crcpmp_pf@meta.data %>%
  dplyr::group_by(manual_l4, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l4), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(unique(seurat_hc_crcpmp_pf@meta.data[,c("manual_l1", "manual_l4")]), by = "manual_l4") %>%
  dplyr::mutate(manual_l1 = factor(manual_l1, levels = manual_l1_order)) %>%
  dplyr::left_join(dacs_pf_l4rl0_crcpmpvhc[[1]]$da, by = c("manual_l4" = "BaselineProp.clusters"))

boxplot_crcpmpvhc_pf_l4rl0_together_ggplotobj <- crcpmpvhc_pf_l4rl0_proportions_df %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l4, -Ncellperc), y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_grid(.~ manual_l1, scales = "free", space="free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_l4rl0_together.pdf", width = 12.5, height = 5)
print(boxplot_crcpmpvhc_pf_l4rl0_together_ggplotobj)
dev.off()

boxplot_crcpmpvhc_pf_l4rl0_separate_ggplotobj <- crcpmpvhc_pf_l4rl0_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l4, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label)),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  ggplot(aes(x = Group, y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%CD45+ cells") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_l4rl0_separate.pdf", width = 13, height = 15)
print(boxplot_crcpmpvhc_pf_l4rl0_separate_ggplotobj)
dev.off()

# Fig boxplot CRC PM+ vs HC PF manual_l3 relative manual_l1

crcpmpvhc_pf_l3rl1_proportions_df <- seurat_hc_crcpmp_pf@meta.data %>%
  dplyr::group_by(manual_l1, manual_l3, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l1, manual_l3), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID, manual_l1) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(dacs_df %>%
                     dplyr::mutate(celltype = rownames(.)), by = c("manual_l3" = "celltype")) %>%
  dplyr::filter(!is.na(P.Value))

boxplot_crcpmpvhc_pf_l3rl1_separate_ggplotobj <- crcpmpvhc_pf_l3rl1_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l3, "\nParent: ", manual_l1, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label)),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  ggplot(aes(x = Group, y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%parent") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_l3rl1_separate.pdf", width = 10, height = 15)
print(boxplot_crcpmpvhc_pf_l3rl1_separate_ggplotobj)
dev.off()

# Fig boxplot CRC PM+ vs HC PF manual_l4 relative manual_l2

crcpmpvhc_pf_l4rl2_proportions_df <- seurat_hc_crcpmp_pf@meta.data %>%
  dplyr::group_by(manual_l2, manual_l4, SampleID, Group, .drop = T) %>% 
  dplyr::summarise(Ncells = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(tidyr::nesting(Group, SampleID), tidyr::nesting(manual_l2, manual_l4), fill = list(Ncells = 0)) %>%
  dplyr::group_by(SampleID, manual_l2) %>%
  dplyr::mutate(Nlrefsample = sum(Ncells),
                Ncellprop = Ncells/Nlrefsample,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Donor = gsub("^(pt[0-9]{2})_.+$", "\\1", SampleID)) %>%
  dplyr::left_join(dacs_df %>%
                     dplyr::mutate(celltype = rownames(.)), by = c("manual_l4" = "celltype")) %>%
  dplyr::filter(!is.na(P.Value))

boxplot_crcpmpvhc_pf_l4rl2_separate_ggplotobj <- crcpmpvhc_pf_l4rl2_proportions_df %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(label = paste0(manual_l4, "\nParent: ", manual_l2, "\np-value = ", formatC(P.Value, format = "e", digits = 3)),
                label = factor(label, levels = unique(label)),
                Group = factor(Group, levels = c("HC", "CRC+"))) %>%
  ggplot(aes(x = Group, y = Ncellperc, col = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  facet_wrap(~label, scales = "free") +
  labs(y = "%parent") +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_crcpmpvhc_pf_l4rl2_separate.pdf", width = 15, height = 17.5)
print(boxplot_crcpmpvhc_pf_l4rl2_separate_ggplotobj)
dev.off()

# Fig dotplot CRC PM+ vs HC PF macrophages manual_l4 degs

crcpmpvhc_pf_macrophages_l4_top50 <- unique(unlist(lapply(degs_pf_l4_crcpmpvhc[macrophage_l4], function(macrophage){
  macrophage$degs[1:50,"gene"]
})))

macrophage_l4_top50_gene_df <- do.call(rbind, lapply(macrophage_l4, function(macrophage){
  degs_pf_l4_crcpmpvhc[[macrophage]]$degs %>%
    data.frame() %>%
    dplyr::filter(gene %in% crcpmpvhc_pf_macrophages_l4_top50) %>%
    dplyr::select(stat, pvalue, padj, gene) %>%
    dplyr::mutate(macrophages = macrophage)
}))


dotplot_crcpmpvhc_pf_l4_degs_top5_ggplotobj <- macrophage_l4_top5_gene_df %>%
  dplyr::mutate(direction = ifelse(stat<0, "HC", "CRC+"),
                significance = ifelse(pvalue<0.05, "Significant", "NS"),
                gene = factor(gene, levels = crcpmpvhc_pf_macrophages_l4_top5),
                macrophages = factor(macrophages, levels = macrophage_l4)) %>%
  ggplot(aes(x = macrophages, y = gene)) +
  geom_point(aes(size = -log10(pvalue), col = direction, alpha = significance)) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("dotplot_crcpmpvhc_pf_l4_degs_top5.pdf", width = 4, height = 8)
print(dotplot_crcpmpvhc_pf_l4_degs_top5_ggplotobj)
dev.off()

# Fig dotplot CRC PM+ vs HC PF macrophages manual_l4 immunosuppressive genes

macrophages_markers <- read.csv("config/genes_of_interest/macrophages.txt", header = F)[,1]
macrophage_set1 <- unique(read.csv("config/genes_of_interest/macrophages_set1.txt", header = F)[,1])
macrophage_set2 <- unique(read.csv("config/genes_of_interest/macrophages_set2.txt", header = F)[,1])
macrophage_set3 <- unique(read.csv("config/genes_of_interest/macrophages_set3.txt", header = F)[,1])
macrophage_set4 <- unique(read.csv("config/genes_of_interest/macrophages_set4.txt", header = F)[,1])

dotplot_crcpmpvhc_pf_l4_goi_ggplotobj <- do.call(rbind, lapply(macrophage_l4, function(macrophage){
  degs_pf_l4_crcpmpvhc[[macrophage]]$degs %>%
    data.frame() %>%
    dplyr::filter(gene %in% macrophage_set4) %>%
    dplyr::select(stat, pvalue, padj, gene) %>%
    dplyr::mutate(macrophages = macrophage)
})) %>%
  dplyr::mutate(direction = ifelse(stat<0, "HC", "CRC+"),
                significance = ifelse(pvalue<0.05, "Significant", "NS"),
                gene = factor(gene, levels = macrophage_set4),
                macrophages = factor(macrophages, levels = macrophage_l4)) %>%
  ggplot(aes(x = macrophages, y = gene)) +
  geom_point(aes(size = -log10(pvalue), col = direction, alpha = significance)) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("dotplot_crcpmpvhc_pf_l4_set1.pdf", width = 4, height = 27.5)
print(dotplot_crcpmpvhc_pf_l4_goi_ggplotobj)
dev.off()

pdf("dotplot_crcpmpvhc_pf_l4_set2.pdf", width = 4, height = 10)
print(dotplot_crcpmpvhc_pf_l4_goi_ggplotobj)
dev.off()

pdf("dotplot_crcpmpvhc_pf_l4_set3.pdf", width = 4, height = 7.5)
print(dotplot_crcpmpvhc_pf_l4_goi_ggplotobj)
dev.off()

pdf("dotplot_crcpmpvhc_pf_l4_set4.pdf", width = 4, height = 12.5)
print(dotplot_crcpmpvhc_pf_l4_goi_ggplotobj)
dev.off()

# Fig heatmap CRC PM+ vs HC PF macrophages manual_l4 fgsea

crcpmpvhc_pf_macrophages_l4_fgsea_sig <- unique(unlist(lapply(fgsea_pf_l4_crcpmpvhc[macrophage_l4], function(macrophage){
  macrophage %>%
    data.frame() %>%
    dplyr::select(1:7) %>%
    dplyr::filter(padj<0.05) %>%
    dplyr::pull(pathway)
})))

crcpmpvhc_pf_macrophages_l4_fgsea_sig_df <- do.call(rbind, lapply(macrophage_l4, function(macrophage){
  fgsea_pf_l4_crcpmpvhc[[macrophage]] %>%
    data.frame() %>%
    dplyr::filter(pathway %in% crcpmpvhc_pf_macrophages_l4_fgsea_sig) %>%
    dplyr::select(pathway, NES, pval, padj) %>%
    dplyr::mutate(macrophages = macrophage)
}))

sig_pwoi <- unique(c("KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_PATHWAYS_IN_CANCER","KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_DNA_REPLICATION","KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_COLORECTAL_CANCER","KEGG_HEDGEHOG_SIGNALING_PATHWAY","KEGG_HEMATOPOIETIC_CELL_LINEAGE","KEGG_REGULATION_OF_ACTIN_CYTOSKELETON","KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_CELL_ADHESION_MOLECULES_CAMS","KEGG_HEMATOPOIETIC_CELL_LINEAGE","KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_TYPE_I_DIABETES_MELLITUS","KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY","KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION","KEGG_TYPE_II_DIABETES_MELLITUS","KEGG_PATHWAYS_IN_CANCER","KEGG_OXIDATIVE_PHOSPHORYLATION","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_PATHWAYS_IN_CANCER","KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_TYPE_II_DIABETES_MELLITUS","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_DNA_REPLICATION","KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_COLORECTAL_CANCER","KEGG_HEDGEHOG_SIGNALING_PATHWAY","KEGG_HEMATOPOIETIC_CELL_LINEAGE","KEGG_REGULATION_OF_ACTIN_CYTOSKELETON"))

dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig_ggplotobj <- crcpmpvhc_pf_macrophages_l4_fgsea_sig_df %>%
  dplyr::mutate(direction = ifelse(NES<0, "HC", "CRC+"),
                significance = ifelse(pval<0.05, "Significant", "NS"),
                pathway = factor(pathway, levels = unique(pathway)),
                macrophages = factor(macrophages, levels = macrophage_l4)) %>%
  ggplot(aes(x = macrophages, y = pathway)) +
  geom_point(aes(size = -log10(pval), col = direction, alpha = significance)) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig.pdf", width = 7.5, height = 15)
print(dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig_ggplotobj)
dev.off()

dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig_pwoi_ggplotobj <- crcpmpvhc_pf_macrophages_l4_fgsea_sig_df %>%
  dplyr::filter(pathway %in% sig_pwoi) %>%
  dplyr::mutate(direction = ifelse(NES<0, "HC", "CRC+"),
                significance = ifelse(pval<0.05, "Significant", "NS"),
                pathway = factor(pathway, levels = unique(pathway)),
                macrophages = factor(macrophages, levels = macrophage_l4)) %>%
  ggplot(aes(x = macrophages, y = pathway)) +
  geom_point(aes(size = -log10(pval), col = direction, alpha = significance)) +
  theme_bw() +
  scale_color_manual(values = group_colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig_pwoi.pdf", width = 6, height = 7.5)
print(dotplot_crcpmpvhc_pf_macrophages_l4_fgsea_sig_pwoi_ggplotobj)
dev.off()

