#!/usr/bin/env Rscript
# This script will create a heatmap of the CD4T cells at the level of celltype l3.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
markers_path <- args[2]
heatmap_svg_path <- args[3]

seuratObject <- readRDS(seurat_rds_path)

# Select the CD4T

cd4t_seuratObject <- seuratObject[,seuratObject@meta.data$manual_l2 == "CD4T"]

# Prepare the markers for CD4T 

cd4t_markers <- readxl::read_excel(markers_path) %>%
  dplyr::filter(celltype_l2 == "CD4T")

cd4t_marker_genes <- paste(paste0(cd4t_markers$Gene_positive, ";", cd4t_markers$Gene_mid, ";", cd4t_markers$Gene_negative), collapse=";")
cd4t_marker_genes <- gsub(" ", "", cd4t_marker_genes)
cd4t_marker_genes <- unique(unlist(strsplit(cd4t_marker_genes, ";")))
cd4t_marker_genes <- cd4t_marker_genes[!cd4t_marker_genes %in% c("[cellcyclegenes]", "NA")]

# Create figure

cd4t_marker_expr <- GetAssayData(cd4t_seuratObject, assay = "RNA")[which(rownames(cd4t_seuratObject) %in% cd4t_marker_genes),]

cd4t_genecelltype_median <- data.frame(GeneID = rownames(cd4t_marker_expr), cd4t_marker_expr) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(cd4t_seuratObject@meta.data), 
                               Donor = cd4t_seuratObject@meta.data$Donor,
                               Celltype = cd4t_seuratObject@meta.data$manual_l3), 
                       by = "CB") %>%
  dplyr::group_by(Celltype, GeneID) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100)

#hclust
cd4t_genecelltype_median_wide <- cd4t_genecelltype_median %>% 
  tidyr::pivot_wider(!Percentage, names_from = Celltype, values_from = Median) %>%
  data.frame()
rownames(cd4t_genecelltype_median_wide) <- cd4t_genecelltype_median_wide$GeneID
cd4t_genecelltype_median_wide <- cd4t_genecelltype_median_wide[,-1]
colnames(cd4t_genecelltype_median_wide) <- gsub("\\.", " ", colnames(cd4t_genecelltype_median_wide))

## Genes
hclust_cd4t_markers <- hclust(dist(cd4t_genecelltype_median_wide), method = "complete", members = NULL)
## Celltypes
hclust_cd4t_celltypes <- hclust(dist(t(cd4t_genecelltype_median_wide)), method = "complete", members = NULL)

cd4t_genecelltype_median$GeneID <- factor(cd4t_genecelltype_median$GeneID, levels = rownames(cd4t_genecelltype_median_wide)[hclust_cd4t_markers$order])
cd4t_genecelltype_median$Celltype <- factor(cd4t_genecelltype_median$Celltype, levels = colnames(cd4t_genecelltype_median_wide)[hclust_cd4t_celltypes$order])

fig_heatmap_markers <- cd4t_genecelltype_median %>% 
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
        panel.border = element_blank())

svglite::svglite(width = 10, height = 4, file = heatmap_svg_path, bg = "white")
print(fig_heatmap_markers)
dev.off()
