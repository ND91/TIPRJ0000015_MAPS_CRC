#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMC and PF samples from HC colored for tissue.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
hc_umap_pbmc_pf_coltissue_pdf <- args[2]
hc_umap_pbmc_pf_coltissue_splittissue_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)
  
umap_df <- data.frame(CB = colnames(seuratObject),
                      Embeddings(seuratObject[["umap"]]),
                      seuratObject@meta.data)

umap_coltissue_ggplot2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, col = Tissue)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
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

pdf(width = 5, height = 5, file = hc_umap_pbmc_pf_coltissue_pdf)
print(umap_coltissue_ggplot2)
dev.off()

umap_coltissue_splittissue_ggplot2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, col = Tissue)) +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = c(PF = "#224FBD", PBMC = "#F30000")) +
  facet_wrap(~Tissue, ncol = 2) +
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

pdf(width = 10, height = 5, file = hc_umap_pbmc_pf_coltissue_splittissue_pdf)
print(umap_coltissue_splittissue_ggplot2)
dev.off()


sessionInfo()