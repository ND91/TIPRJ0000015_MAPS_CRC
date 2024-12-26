#!/usr/bin/env Rscript

# This script generates a figure addressing reviewer 1, point 8.

require(Seurat)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(ggpubr)
require(ggdendro)

seurat_rds <- "output/subsets/live_singlet_nonproliferating_SeuratObject.Rds"

seuratObject <- readRDS(seurat_rds)

cells_selected <- seuratObject@meta.data %>%
  dplyr::filter(Donor %in% c("pt87", "pt88", "pt89", "pt91", "pt92", "pt31", "pt35", "pt37", "pt71", "pt73", "pt74", "pt76", "pt78"),
                manual_l2 %in% c("Macrophages", "Monocytes")) %>%
  dplyr::pull(CellID) %>%
  unique()

seuratObject <- seuratObject[,which(seuratObject@meta.data$CellID %in% cells_selected)]

seuratObject <- DietSeurat(seuratObject, counts = T, data = T, scale.data = F)
seuratObject <- seuratObject[Matrix::rowSums(seuratObject) != 0, ]
options(future.globals.maxSize= 20000*1024^2)
seuratObject <- SCTransform(seuratObject, conserve.memory = T)
seuratObject <- RunPCA(object = seuratObject, npcs = 100, seed.use = 326126)
seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:42)
seuratObject <- FindClusters(seuratObject, resolution = 1, verbose = FALSE)
seuratObject <- RunUMAP(seuratObject, dims = 1:42, seed.use = 123632)

seuratObject@meta.data$Tissue_group <- factor(paste0(seuratObject@meta.data$Tissue, " ", seuratObject@meta.data$Group), levels = c("PBMC HC", "PBMC CRC+", "PF HC", "PF CRC+", "TX CRC+"))

figR1P8A <- DimPlot(seuratObject, group.by = "Tissue_group")
figR1P8B <- DimPlot(seuratObject, group.by = "manual_l2", label = T)
figR1P8C <- DimPlot(seuratObject, group.by = "Tissue", split.by = "Group")

seuratObject@meta.data$Tissue_group_donor_l2 <- paste0(seuratObject@meta.data$Tissue, " ", seuratObject@meta.data$Group, " ", seuratObject@meta.data$Donor, " ", seuratObject@meta.data$manual_l2)

counts_tissue_group <- GetAssayData(seuratObject, assay = "RNA", layer = "counts") %*% model.matrix(~0+Tissue_group_donor_l2, data = seuratObject@meta.data)
colnames(counts_tissue_group) <- gsub("Tissue_group_donor_l2", "", colnames(counts_tissue_group))

counts_tissue_group <- counts_tissue_group[which(Matrix::rowSums(counts_tissue_group) != 0),]

rld_counts_tissue_group <- DESeq2::rlog(as.matrix(counts_tissue_group))

dist_counts_tissue_group <- dist(t(rld_counts_tissue_group), method = 'euclidean')

figR1P8D <- as.ggplot(pheatmap(as.matrix(dist_counts_tissue_group), display_numbers = T, number_format = `%.0f`))

figR1P8AB <- ggarrange(figR1P8A, figR1P8C, nrow = 1, ncol = 2, widths = c(1.25, 2))
figR1P8CD <- ggarrange(figR1P8B, figR1P8D, nrow = 1, ncol = 2, widths = c(1.25, 2))

figR1P8 <- ggarrange(figR1P8AB, figR1P8CD, nrow = 2, ncol = 1)

ggsave(filename = "figR1P8.pdf", width = 18, height = 12, units = "in")
