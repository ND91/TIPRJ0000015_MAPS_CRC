#!/usr/bin/env Rscript
# This script will be used to quantify and characterize the macrophages observed in PBMCs.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(svglite))

args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 5) {
#   stop(paste0("Script needs 5 arguments. Current input is:", args))
# }

seurat_myeloid_rds_path <- args[1]
seurat_macrophages_rds_path <- args[2]

myeloid_seuratObject <- readRDS(seurat_myeloid_rds_path)

# Figure 1: tSNE colored by manual_l2, split by tissue

myeloid_seuratObject_metadata_tsne <- data.frame(Embeddings(myeloid_seuratObject[["tsne"]]),
                                                 myeloid_seuratObject@meta.data) %>%
  dplyr::mutate(manual_l2_number = as.numeric(as.factor(manual_l2)),
                manual_l2_w_number = paste0(as.numeric(as.factor(manual_l2)), ". ", manual_l2),
                manual_l3_number = as.numeric(as.factor(manual_l3)),
                manual_l3_w_number = paste0(as.numeric(as.factor(manual_l3)), ". ", manual_l3),
                manual_l4_number = as.numeric(as.factor(manual_l4)),
                manual_l4_w_number = paste0(as.numeric(as.factor(manual_l4)), ". ", manual_l4),
                Group1 = factor(ifelse(Group == "HC", "HC", "CRC"), levels = c("CRC", "HC")),
                Group2 = factor(Group, levels = c("HC", "CRC-", "CRC+")))

fig1 <- myeloid_seuratObject_metadata_tsne %>%
  dplyr::filter(Group != "HC") %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2)) +
  geom_point_rast(show.legend = T, size = 2, col = "black") +
  geom_point_rast(show.legend = T, size = 1, aes(col = manual_l2_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = myeloid_seuratObject_metadata_tsne %>%
                     dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                     summarize(x = median(x = tSNE_1),
                               y = median(x = tSNE_2)),
                   mapping = aes(label = manual_l2_number, x = x, y = y, col = manual_l2_w_number),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(subtitle = "PF and PBMCs from CRC patients") +
  facet_grid(.~Tissue) +
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
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"))

svglite(width = 10, height = 5, file = fig1_svg_path, bg = "white")
print(fig1)
dev.off()

# Figure 2: tSNE superposed by expression of myeloid markers CD163, MARCO, CD68, CD14, FCGR3A, CLEC9A

myeloid_markers <- c("CD14", "FCGR3A", "CD163", "MARCO", "CD68", "CD1C")

myeloid_seuratObject_marker_tsne <- data.frame(Embeddings(myeloid_seuratObject[["tsne"]]),
                                               myeloid_seuratObject@meta.data,
                                               t(as.matrix(GetAssayData(myeloid_seuratObject)[myeloid_markers,]))) %>%
  dplyr::mutate(manual_l2_number = as.numeric(as.factor(manual_l2)),
                manual_l2_w_number = paste0(as.numeric(as.factor(manual_l2)), ". ", manual_l2),
                manual_l3_number = as.numeric(as.factor(manual_l3)),
                manual_l3_w_number = paste0(as.numeric(as.factor(manual_l3)), ". ", manual_l3),
                manual_l4_number = as.numeric(as.factor(manual_l4)),
                manual_l4_w_number = paste0(as.numeric(as.factor(manual_l4)), ". ", manual_l4),
                Group1 = factor(ifelse(Group == "HC", "HC", "CRC"), levels = c("CRC", "HC")),
                Group2 = factor(Group, levels = c("HC", "CRC-", "CRC+")))

fig2_list <- lapply(myeloid_markers, function(marker_gene){
  myeloid_seuratObject_marker_tsne %>%
    dplyr::filter(Group != "HC") %>%
    ggplot(aes(x = tSNE_1, y = tSNE_2)) +
    geom_point_rast(show.legend = F, size = 2, col = "black") +
    geom_point_rast(show.legend = T, size = 1, aes_string(col = marker_gene, alpha = marker_gene)) +
    labs(y = "",
         x = "") +
    geom_label_repel(data = myeloid_seuratObject_marker_tsne %>%
                       dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                       summarize(x = median(x = tSNE_1),
                                 y = median(x = tSNE_2)),
                     mapping = aes(label = manual_l2_number, x = x, y = y),
                     alpha = 1, 
                     show.legend = F) +
    #guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(subtitle = marker_gene) +
    #facet_wrap(~Tissue) +
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
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold"))
})

fig2 <- ggarrange(plotlist = fig2_list, nrow = 2, ncol = 3)

svglite(width = 15, height = 10, file = fig2_svg_path, bg = "white")
print(fig2)
dev.off()

# Figure 3: tSNE colored by manual_l2, for PBMC alone from HC and CRC patients. Split by group.

fig3 <- myeloid_seuratObject_metadata_tsne %>%
  dplyr::filter(Tissue == "PBMC") %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2)) +
  geom_point_rast(show.legend = T, size = 2, col = "black") +
  geom_point_rast(show.legend = T, size = 1, aes(col = manual_l2_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = myeloid_seuratObject_metadata_tsne %>%
                     dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                     summarize(x = median(x = tSNE_1),
                               y = median(x = tSNE_2)),
                   mapping = aes(label = manual_l2_number, x = x, y = y, col = manual_l2_w_number),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(subtitle = "PBMCs from HC and CRC patients") +
  facet_wrap(~Group2, nrow = 1, ncol = 3) +
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

svglite(width = 15, height = 5, file = fig3_svg_path, bg = "white")
print(fig3)
dev.off()

# Figure 4: Boxplot overlaid with a jitterplot representing the abundance of the macrophages (l2 relative to l1) per tissue.

macrophages_CRCvHC_da <- myeloid_seuratObject_metadata_tsne %>%
  dplyr::filter(Tissue == "PBMC") %>%
  propeller(sample = .$SampleID, clusters = .$manual_l2, group = .$Group1)

macrophages_CRCpvCRCn_da <- myeloid_seuratObject_metadata_tsne %>%
  dplyr::filter(Tissue == "PBMC",
                Group != "HC") %>%
  dplyr::mutate(Group2 = as.factor(as.character(Group2))) %>%
  propeller(sample = .$SampleID, clusters = .$manual_l2, group = .$Group2)

fig4 <- data.frame(table(myeloid_seuratObject_metadata_tsne$manual_l2, myeloid_seuratObject_metadata_tsne$SampleID)) %>%
  dplyr::rename(manual_l2 = Var1,
                SampleID = Var2) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Total = sum(Freq),
                Proportion = Freq/Total) %>%
  dplyr::left_join(unique(myeloid_seuratObject_metadata_tsne[,c("SampleID", "Group", "Tissue")]), by = "SampleID") %>%
  dplyr::filter(manual_l2 == "Macrophages",
                Tissue == "PBMC") %>%
  dplyr::mutate(Group = factor(Group, levels = c("HC", "CRC-", "CRC+"))) %>%
  ggplot(aes(x = Group, y = Proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point_rast(alpha = 0.5) +
  labs(title = "Proportion macrophages relative to myeloid",
       subtitle = paste0("p-value CRC vs HC = ", formatC(macrophages_CRCvHC_da$P.Value[1], digits = 2, format = "e"), "\np-value CRC+ vs CRC- = ", round(macrophages_CRCpvCRCn_da$P.Value[3], 3))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

svglite(width = 5, height = 5, file = fig4_svg_path, bg = "white")
print(fig4)
dev.off()

# Figure 5: tSNE of the macrophages, colored by macrophage markers

macrophage_markers <- c("CD163", "C1QC", "FN1", "ACTB", "SPP1", "VCAN", "CD163", "GATA6")

macrophage_seuratObject_marker_tsne <- data.frame(Embeddings(macrophages_seuratObject[["tsne"]]),
                                                  macrophages_seuratObject@meta.data,
                                               t(as.matrix(GetAssayData(macrophages_seuratObject)[macrophage_markers,]))) %>%
  dplyr::mutate(manual_l2_number = as.numeric(as.factor(manual_l2)),
                manual_l2_w_number = paste0(as.numeric(as.factor(manual_l2)), ". ", manual_l2),
                manual_l3_number = as.numeric(as.factor(manual_l3)),
                manual_l3_w_number = paste0(as.numeric(as.factor(manual_l3)), ". ", manual_l3),
                manual_l4_number = as.numeric(as.factor(manual_l4)),
                manual_l4_w_number = paste0(as.numeric(as.factor(manual_l4)), ". ", manual_l4),
                Group1 = factor(ifelse(Group == "HC", "HC", "CRC"), levels = c("CRC", "HC")),
                Group2 = factor(Group, levels = c("HC", "CRC-", "CRC+")))

fig5_list <- lapply(macrophage_markers, function(marker_gene){
  macrophage_seuratObject_marker_tsne %>%
    dplyr::filter(Group != "HC") %>%
    ggplot(aes(x = tSNE_1, y = tSNE_2)) +
    geom_point_rast(show.legend = F, size = 2, col = "black") +
    geom_point_rast(show.legend = T, size = 1, aes_string(col = marker_gene, alpha = marker_gene)) +
    labs(y = "",
         x = "") +
    geom_label_repel(data = myeloid_seuratObject_marker_tsne %>%
                       dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                       summarize(x = median(x = tSNE_1),
                                 y = median(x = tSNE_2)),
                     mapping = aes(label = manual_l2_number, x = x, y = y),
                     alpha = 1, 
                     show.legend = F) +
    #guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(subtitle = marker_gene) +
    #facet_wrap(~Tissue) +
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
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold"))
})

fig5 <- ggarrange(plotlist = fig2_list, nrow = 2, ncol = 4)

svglite(width = 20, height = 10, file = fig5_svg_path, bg = "white")
print(fig5)
dev.off()

# Figure 6: tSNE of the macrophages colored by manual_l4, split by tissue

macrophages_seuratObject_metadata_tsne <- data.frame(Embeddings(macrophages_seuratObject[["tsne"]]),
                                                     macrophages_seuratObject@meta.data) %>%
  dplyr::mutate(manual_l2_number = as.numeric(as.factor(manual_l2)),
                manual_l2_w_number = paste0(as.numeric(as.factor(manual_l2)), ". ", manual_l2),
                manual_l3_number = as.numeric(as.factor(manual_l3)),
                manual_l3_w_number = paste0(as.numeric(as.factor(manual_l3)), ". ", manual_l3),
                manual_l4_number = as.numeric(as.factor(manual_l4)),
                manual_l4_w_number = paste0(as.numeric(as.factor(manual_l4)), ". ", manual_l4),
                Group1 = factor(ifelse(Group == "HC", "HC", "CRC"), levels = c("CRC", "HC")),
                Group2 = factor(Group, levels = c("HC", "CRC-", "CRC+")))

fig6 <- macrophages_seuratObject_metadata_tsne %>%
  #dplyr::filter(Group != "HC") %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2)) +
  geom_point_rast(show.legend = T, size = 2, col = "black") +
  geom_point_rast(show.legend = T, size = 1, aes(col = manual_l4_w_number)) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = macrophages_seuratObject_metadata_tsne %>%
                     dplyr::group_by(manual_l4_number, manual_l4_w_number) %>%
                     summarize(x = median(x = tSNE_1),
                               y = median(x = tSNE_2)),
                   mapping = aes(label = manual_l4_number, x = x, y = y, col = manual_l4_w_number),
                   alpha = 1, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(subtitle = "Macrophages") +
  facet_grid(.~Tissue) +
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
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"))

svglite(width = 10, height = 5, file = fig6_svg_path, bg = "white")
print(fig6)
dev.off()

# Figure 7: Boxplot

macrophages_subsets_CRCpvCRCn_da <- macrophages_seuratObject_metadata_tsne %>%
  dplyr::filter(Tissue == "PBMC",
                Group != "HC") %>%
  dplyr::mutate(Group2 = as.factor(as.character(Group2))) %>%
  propeller(sample = .$SampleID, clusters = .$manual_l4, group = .$Group2)

fig7 <- data.frame(table(macrophages_seuratObject_metadata_tsne$manual_l4, macrophages_seuratObject_metadata_tsne$SampleID)) %>%
  dplyr::rename(manual_l4 = Var1,
                SampleID = Var2) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Total = sum(Freq),
                Proportion = Freq/Total) %>%
  dplyr::left_join(unique(macrophages_seuratObject_metadata_tsne[,c("SampleID", "Group", "Tissue")]), by = "SampleID") %>%
  dplyr::mutate(Group = factor(Group, levels = c("CRC-", "CRC+"))) %>%
  dplyr::filter(Tissue == "PBMC") %>%
  dplyr::left_join(macrophages_subsets_CRCpvCRCn_da, by = c("manual_l4" = "BaselineProp.clusters")) %>%
  dplyr::mutate(label = paste0(manual_l4, "\np-value = ", round(P.Value, 3))) %>%
  ggplot(aes(x = Group, y = Proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point_rast(alpha = 0.5) +
  labs(title = "Proportion PBMC macrophages subsets relative to all macrophages") +
  facet_wrap(~label, nrow = 2, ncol = 4, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

svglite(width = 10, height = 5, file = fig7_svg_path, bg = "white")
print(fig7)
dev.off()
