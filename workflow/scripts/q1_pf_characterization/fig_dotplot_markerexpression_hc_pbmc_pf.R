#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l1 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste0("Script needs 7 arguments. Current input is:", args))
}

seurat_rds <- args[1]
markers_xlsx <- args[2]
celltype_level <- args[3]
dotplot_pdf <- args[4]
dotplot_splittissue_pdf <- args[5]
dotplot_pf_gex_pdf <- args[6]
dotplot_pf_pex_pdf <- args[7]

seuratObject <- readRDS(seurat_rds)

marker_genes <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% unique(seuratObject@meta.data[,celltype_level]),
                modality == "gene")

marker_proteins <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% unique(seuratObject@meta.data[,celltype_level]),
                modality == "protein") %>%
  dplyr::mutate(canonical_marker = gsub("^Hu\\.", "", canonical_marker))

marker_gex <- GetAssayData(seuratObject, assay = "RNA")[which(rownames(GetAssayData(seuratObject, assay = "RNA")) %in% unique(marker_genes$canonical_marker)), ]

# Together

marker_gex_median <- data.frame(GeneID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data[,celltype_level]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, unique(marker_genes$celltype)),
                GeneID = factor(GeneID, unique(marker_genes$canonical_marker)))

dotplot <- marker_gex_median %>% 
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = GeneID, y = Celltype, col = Median)) +
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

if(celltype_level == "manual_l1"){
  pdf(width = 6, height = 2.9, file = dotplot_pdf)
  print(dotplot)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 8.5, height = 5, file = dotplot_pdf)
  print(dotplot)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 12, height = 8, file = dotplot_pdf)
  print(dotplot)
  dev.off()
}

# Split per tissue

marker_gex_median_splittissue <- data.frame(GeneID = rownames(marker_gex), marker_gex) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data),
                               Tissue = seuratObject@meta.data[,"Tissue"],
                               Celltype = seuratObject@meta.data[,celltype_level]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID),
                Tissue = factor(Tissue)) %>%
  dplyr::group_by(Celltype, GeneID, Tissue, .drop = T) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, unique(marker_genes$celltype)),
                GeneID = factor(GeneID, unique(marker_genes$canonical_marker)),
                Celltype_Tissue = paste0(Celltype, " (", Tissue, ")"))

marker_gex_median_splittissue_order <- data.frame(expand.grid(levels(marker_gex_median_splittissue$Celltype), unique(marker_gex_median_splittissue$Tissue))) %>%
  dplyr::rename(Celltype = Var1,
                Tissue = Var2) %>%
  dplyr::mutate(Celltype_Tissue = paste0(Celltype, " (", Tissue, ")")) %>%
  dplyr::arrange(Celltype)

marker_gex_median_splittissue$Celltype_Tissue <- factor(marker_gex_median_splittissue$Celltype_Tissue, levels = marker_gex_median_splittissue_order$Celltype_Tissue)

dotplot_splittissue <- marker_gex_median_splittissue %>% 
  ggplot(aes(x = Celltype_Tissue, y = GeneID, col = Median)) +
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

if(celltype_level == "manual_l1"){
  pdf(width = 5, height = 9, file = dotplot_splittissue_pdf)
  print(dotplot_splittissue)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 7.5, height = 10, file = dotplot_splittissue_pdf)
  print(dotplot_splittissue)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 12, height = 13, file = dotplot_splittissue_pdf)
  print(dotplot_splittissue)
  dev.off()
}

# PF gex

marker_gex_median_pf <- marker_gex_median_splittissue %>%
  dplyr::filter(Tissue == "PF")

dotplot_pf_gex <- marker_gex_median_pf %>% 
  ggplot(aes(x = Celltype, y = GeneID, col = Median)) +
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

if(celltype_level == "manual_l1"){
  pdf(width = 1.75, height = 5, file = dotplot_pf_gex_pdf)
  print(dotplot_pf_gex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 4.5, height = 10, file = dotplot_pf_gex_pdf)
  print(dotplot_pf_gex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 7.5, height = 15, file = dotplot_pf_gex_pdf)
  print(dotplot_pf_gex)
  dev.off()
}

# PF pex

seuratObject_pf <- seuratObject[,seuratObject@meta.data$Tissue == "PF"]
DefaultAssay(seuratObject_pf) <- "CITE"

marker_pex <- GetAssayData(seuratObject_pf, assay = "CITE")[which(rownames(GetAssayData(seuratObject_pf, assay = "CITE")) %in% paste0("Hu.", unique(marker_proteins$canonical_marker))), ]

marker_pex_median <- data.frame(ProteinID = gsub("^Hu\\.", "", rownames(marker_pex)), marker_pex) %>%
  tidyr::pivot_longer(-ProteinID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject_pf@meta.data), 
                               Celltype = seuratObject_pf@meta.data[,celltype_level]), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                ProteinID = factor(ProteinID)) %>%
  dplyr::group_by(Celltype, ProteinID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, unique(marker_proteins$celltype)),
                ProteinID = factor(ProteinID, unique(marker_proteins$canonical_marker)))

dotplot_pf_pex <- marker_pex_median %>% 
  #dplyr::filter(Percentage != 0) %>%
  ggplot(aes(x = ProteinID, y = Celltype, col = Median)) +
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

if(celltype_level == "manual_l1"){
  pdf(width = 6, height = 3.25, file = dotplot_pf_pex_pdf)
  print(dotplot_pf_pex)
  dev.off()
} else if(celltype_level == "manual_l2"){
  pdf(width = 6.5, height = 5, file = dotplot_pf_pex_pdf)
  print(dotplot_pf_pex)
  dev.off()
} else if(celltype_level == "manual_l3"){
  pdf(width = 7.5, height = 8, file = dotplot_pf_pex_pdf)
  print(dotplot_pf_pex)
  dev.off()
}

sessionInfo()
