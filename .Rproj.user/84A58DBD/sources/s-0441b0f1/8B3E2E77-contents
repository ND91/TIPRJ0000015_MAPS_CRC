#!/usr/bin/env Rscript
# This script will generate the following plots from the PBMC, PF, and TX tissue samples obtained from donors that provided all.
# 1. A tSNE plot of all cells together colored by tissue.
# 2. A tSNE plot of all cells together facetted by tissue, colored by manual_l2.
# 3. A tSNE plot of all cells facetted by manual_l1, colored by manual_l3.

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(svglite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
plot_df_csv_path <- args[2]
fig1_svg_path <- args[3]
fig2_svg_path <- args[4]
fig3_svg_path <- args[5]

seuratObject <- readRDS(seurat_rds_path)

seuratObject_metadata_tsne <- data.frame(Embeddings(seuratObject[["tsne"]]),
                                         seuratObject@meta.data) %>%
  dplyr::mutate(manual_l2_number = as.numeric(as.factor(manual_l2)),
                manual_l2_w_number = paste0(as.numeric(as.factor(manual_l2)), ". ", manual_l2),
                manual_l3_number = as.numeric(as.factor(manual_l3)),
                manual_l3_w_number = paste0(as.numeric(as.factor(manual_l3)), ". ", manual_l3))

write.csv(seuratObject_metadata_tsne, plot_df_csv_path)

# 1. A tSNE plot of all cells together colored by tissue.

fig1 <- ggplot(seuratObject_metadata_tsne, aes(x = tSNE_1, y = tSNE_2, col = Tissue)) +
  geom_point_rast(show.legend = T, size = 0.1) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
        )

svglite(width = 7.5, height = 7.5, file = fig1_svg_path, bg = "white")
print(fig1)
dev.off()

# 2. A tSNE plot of all cells together facetted by tissue, colored by manual_l2.

fig2 <- ggplot(seuratObject_metadata_tsne, aes(x = tSNE_1, y = tSNE_2, col = manual_l2_w_number)) +
  geom_point_rast(show.legend = T, size = 0.1) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = seuratObject_metadata_tsne %>%
                     dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                     summarize(x = median(x = tSNE_1),
                               y = median(x = tSNE_2)),
                   mapping = aes(label = manual_l2_number, x = x, y = y),
                   alpha = 0.9, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~Tissue, nrow = 1, ncol = 3) +
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

svglite(width = 22.5, height = 7.5, file = fig2_svg_path, bg = "white")
print(fig2)
dev.off()

# 3. A tSNE plot per manual_l1, colored by manual_l3.

fig3 <- seuratObject_metadata_tsne %>%
  dplyr::filter(manual_l1 %in% c("T", "B", "NK/ILC", "Myeloid")) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, col = manual_l2_w_number)) +
  geom_point_rast(show.legend = T, size = 0.1) +
  labs(y = "",
       x = "") +
  geom_label_repel(data = seuratObject_metadata_tsne %>%
                     dplyr::group_by(manual_l2_number, manual_l2_w_number) %>%
                     summarize(x = median(x = tSNE_1),
                               y = median(x = tSNE_2)),
                   mapping = aes(label = manual_l2_number, x = x, y = y),
                   alpha = 0.9, 
                   show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~manual_l1, nrow = 1, ncol = 4) +
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

svglite(width = 4*7.5, height = 7.5, file = fig3_svg_path, bg = "white")
print(fig3)
dev.off()