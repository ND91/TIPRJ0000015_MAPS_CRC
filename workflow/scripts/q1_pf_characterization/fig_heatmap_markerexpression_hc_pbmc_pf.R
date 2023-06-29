#!/usr/bin/env Rscript

# The goal of this script is to create a heatmap of the markergenes (and proteins in the case of PF) of all PBMC and PF samples from HC.

library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0("Script needs 8 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_xlsx <- args[2]
celltype_level <- args[3]
heatmap_pdf <- args[4]
heatmap_splittissue_pdf <- args[5]
heatmap_pf_gex_pdf <- args[6]
heatmap_pf_pex_pdf <- args[7]
heatmap_pf_gex_pex_pdf <- args[8]

seuratObject <- readRDS(seurat_rds)

marker_genes <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% unique(seuratObject@meta.data[,celltype_level]),
                modality == "gene")

marker_proteins <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% unique(seuratObject@meta.data[,celltype_level]),
                modality == "protein")

# Together

Idents(seuratObject) <- celltype_level
seuratObject_avexpr <- AverageExpression(seuratObject, return.seurat = T)

ol_markers <- unique(marker_genes$canonical_marker)[unique(marker_genes$canonical_marker) %in% rownames(GetAssayData(seuratObject_avexpr, assay = "RNA"))]

marker_avegexpr <- GetAssayData(seuratObject_avexpr, assay = "RNA")[ol_markers, ]

heatmap_annotation <- data.frame(Celltype = rownames(seuratObject_avexpr@meta.data),
                                 row.names = colnames(seuratObject_avexpr)) %>% 
  dplyr::filter(Celltype %in% marker_genes$celltype) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = unique(marker_genes$celltype))) %>%
  dplyr::arrange(Celltype)

marker_avegexpr <- marker_avegexpr[ol_markers, rownames(heatmap_annotation)]
heatmap_gex <- ComplexHeatmap::pheatmap(marker_avegexpr, 
                                        color = PurpleAndYellow(1000),
                                        cluster_rows = F, 
                                        cluster_cols = F, 
                                        scale = "row",
                                        name = "gex")

if(celltype_level == "manual_l1"){
  pdf(width = 2, height = 3, file = heatmap_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 3.5, height = 7.5, file = heatmap_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 6, height = 9, file = heatmap_pdf)
  print(heatmap_gex)
  dev.off()
}

# Split per tissue

seuratObject@meta.data[,celltype_level] <- factor(seuratObject@meta.data[,celltype_level], levels = unique(marker_genes$celltype))

seuratObject@meta.data$Celltype_Tissue = paste0(seuratObject@meta.data[,celltype_level], " (", seuratObject@meta.data$Tissue, ")")
Idents(seuratObject) <- "Celltype_Tissue"
seuratObject_avexpr_splittissue <- AverageExpression(seuratObject, return.seurat = T)

ol_markers_splittissue <- intersect(marker_genes$canonical_marker, rownames(GetAssayData(seuratObject_avexpr_splittissue, assay = "RNA")))

marker_avegexpr_splittissue <- GetAssayData(seuratObject_avexpr_splittissue, assay = "RNA")[ol_markers_splittissue, ]

marker_expr_median_splittissue_order <- data.frame(expand.grid(levels(seuratObject@meta.data[,celltype_level]), unique(seuratObject@meta.data$Tissue))) %>%
  dplyr::rename(Celltype = Var1,
                Tissue = Var2) %>%
  dplyr::mutate(Celltype_Tissue = paste0(Celltype, " (", Tissue, ")")) %>%
  dplyr::arrange(Celltype) %>%
  dplyr::filter(Celltype_Tissue %in% colnames(marker_avegexpr_splittissue))

marker_avegexpr_splittissue <- marker_avegexpr_splittissue[ol_markers_splittissue, marker_expr_median_splittissue_order$Celltype_Tissue]

heatmap_gex <- ComplexHeatmap::pheatmap(marker_avegexpr_splittissue, 
                                        color = PurpleAndYellow(1000),
                                        cluster_rows = F, 
                                        cluster_cols = F, 
                                        scale = "row",
                                        name = "gex")

if(celltype_level == "manual_l1"){
  pdf(width = 2.75, height = 4, file = heatmap_splittissue_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 6, height = 7.5, file = heatmap_splittissue_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 9, height = 9, file = heatmap_splittissue_pdf)
  print(heatmap_gex)
  dev.off()
}

# PF gex

marker_expr_median_pf_gex_order <- marker_expr_median_splittissue_order %>%
  dplyr::filter(Tissue == "PF")

marker_avegexpr_pf <- marker_avegexpr_splittissue[ol_markers_splittissue, marker_expr_median_pf_gex_order$Celltype_Tissue]
colnames(marker_avegexpr_pf) <- marker_expr_median_pf_gex_order$Celltype

heatmap_gex <- ComplexHeatmap::pheatmap(marker_avegexpr_pf, 
                                        color = PurpleAndYellow(1000),
                                        cluster_rows = F, 
                                        cluster_cols = F, 
                                        scale = "row",
                                        name = "gex")

if(celltype_level == "manual_l1"){
  pdf(width = 2.75, height = 4, file = heatmap_pf_gex_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 6, height = 7.5, file = heatmap_pf_gex_pdf)
  print(heatmap_gex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 9, height = 9, file = heatmap_pf_gex_pdf)
  print(heatmap_gex)
  dev.off()
}

# PF pex

seuratObject_pf <- seuratObject[,seuratObject@meta.data$Tissue == "PF"]
DefaultAssay(seuratObject_pf) <- "CITE"
seuratObject_pf <- NormalizeData(seuratObject_pf)

seuratObject_pf@meta.data[,celltype_level] <- factor(seuratObject_pf@meta.data[,celltype_level], levels = unique(marker_proteins$celltype))

Idents(seuratObject_pf) <- celltype_level
seuratObject_pf_avexpr <- AverageExpression(seuratObject_pf, return.seurat = T)

ol_markers_pex <- intersect(marker_proteins$canonical_marker, rownames(GetAssayData(seuratObject_pf_avexpr, assay = "CITE")))

marker_avepexpr_pf <- GetAssayData(seuratObject_pf_avexpr, assay = "CITE")[ol_markers_pex, ]

heatmap_pex <- ComplexHeatmap::pheatmap(marker_avepexpr_pf, 
                                        #color = colorRampPalette(c("red", "yellow", "blue"))(1000),
                                        color = PurpleAndYellow(1000),
                                        labels_row = gsub("^Hu\\.", "", rownames(marker_avepexpr_pf)),
                                        cluster_rows = F, 
                                        cluster_cols = F, 
                                        scale = "row",
                                        name = "pex")

if(celltype_level == "manual_l1"){
  pdf(width = 2.75, height = 4, file = heatmap_pf_pex_pdf)
  print(heatmap_pex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 6, height = 7.5, file = heatmap_pf_pex_pdf)
  print(heatmap_pex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 9, height = 9, file = heatmap_pf_pex_pdf)
  print(heatmap_pex)
  dev.off()
}

# PF gex + pex

heatmap_gex_pex <- heatmap_gex %v% heatmap_pex

if(celltype_level == "manual_l1"){
  pdf(width = 2.5, height = 8, file = heatmap_pf_gex_pex_pdf)
  print(heatmap_gex_pex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 3.5, height = 8.5, file = heatmap_pf_gex_pex_pdf)
  print(heatmap_gex_pex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 6.5, height = 15, file = heatmap_pf_gex_pex_pdf)
  print(heatmap_gex_pex)
  dev.off()
}


sessionInfo()