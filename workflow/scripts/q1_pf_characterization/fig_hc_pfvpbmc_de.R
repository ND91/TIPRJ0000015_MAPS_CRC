#!/usr/bin/env Rscript

# The goal of this script is to create an arrowplot of all the significant genes identified as differentially expressed between PF and PBMC celltypes.

require(ggplot2)
require(dplyr)
require(ggrastr)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

degs_list_rds <- args[1]
celltype_level <- args[2]
markers_xlsx <- args[3]
arrowplot_top5sig_pdf <- args[4]
arrowplot_resident_pdf <- args[5]

degs_list <- readRDS(degs_list_rds)

celltype_order <- readxl::read_excel(markers_xlsx) %>%
  dplyr::filter(level == celltype_level,
                celltype %in% names(degs_list)) %>%
  dplyr::pull(celltype) %>%
  unique()

# Arrowplot

degs <- lapply(names(degs_list), function(celltype){
  data.frame(degs_list[[celltype]]$degs, celltype = celltype)
})

degs_df <- do.call(rbind, degs) %>%
  dplyr::mutate(global_padj = p.adjust(pvalue, method = "BH"))

# # DEG selection: overlapping DEGs
# 
# degs_sig_df <- degs_df %>%
#   dplyr::filter(padj < 0.05)
# 
# ## I cannot plot all the 10832, some filtering needs to be done. I can however plot genes that are present in 2 or more studies. This will still yield 3246 genes. At 5 or more studies, we are talking about 30 genes, which are a lot more manageable.
# genes_of_interest <- names(which(table(degs_sig_df$gene) >= 5))
# 

# DEG selection: top 5 significant genes per celltype

top5sig_genes_of_interest <- degs_df %>%
  dplyr::mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::group_by(celltype) %>%
  dplyr::arrange(celltype) %>%
  dplyr::slice_max(order_by = rev(pvalue), n = 5) %>%
  dplyr::pull(gene) %>%
  unique()

top5sig_goi_df <- degs_df %>%
  dplyr::filter(gene %in% genes_of_interest) %>%
  dplyr::mutate(celltype = factor(celltype, levels = celltype_order),
                gene = factor(gene, levels = genes_of_interest))

arrowplot_top5sig_ggobj <- goi_df %>% 
  dplyr::mutate(Status = factor(ifelse(stat>0, "Up", "Down"), levels = c("Up", "Down")),
                Significant = ifelse(padj<0.05, "Significant", "NS")) %>%
  ggplot(aes(x = celltype, y = gene)) +
  geom_point(aes(colour = Status,
                 shape = Status,
                 fill = log2FoldChange,
                 size = -log10(pvalue),
                 alpha = Significant)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  scale_color_manual(values = c("coral3","deepskyblue3")) +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(shape=17))) +
  theme(axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

if(celltype_level == "manual_l2"){
  pdf(file = arrowplot_top5sig_pdf, width = 7, height = 12)
} else if(celltype_level == "manual_l3"){
  pdf(file = arrowplot_top5sig_pdf, width = 7, height = 15)
}
print(arrowplot_top5sig_ggobj)
dev.off()

# DEG selection: Residency

#residency_genes_of_interest <- c("CD69", "ITGAE", "ADORA1", "ADORA2A", "ADORA2B", "ADORA3", "CD38", "ENTPD1", "NT5E")
residency_genes_of_interest <- c("CD69", "ITGAE")

residency_goi_df <- degs_df %>%
  dplyr::filter(gene %in% residency_genes_of_interest) %>%
  dplyr::mutate(celltype = factor(celltype, levels = celltype_order),
                gene = factor(gene, levels = residency_genes_of_interest))

arrowplot_residency_ggobj <- residency_goi_df %>% 
  dplyr::mutate(Status = factor(ifelse(stat>0, "Up", "Down"), levels = c("Up", "Down")),
                Significant = ifelse(padj<0.05, "Significant", "NS")) %>%
  ggplot(aes(x = celltype, y = gene)) +
  geom_point(aes(colour = Status,
                 shape = Status,
                 fill = log2FoldChange,
                 size = -log10(pvalue),
                 alpha = Significant)) +
  scale_alpha_discrete(range = c(0.25, 1)) +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  scale_color_manual(values = c("coral3","deepskyblue3")) +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(shape=17))) +
  theme(axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

if(celltype_level == "manual_l2"){
  pdf(file = arrowplot_resident_pdf, width = 5, height = 2)
} else if(celltype_level == "manual_l3"){
  pdf(file = arrowplot_resident_pdf, width = 7, height = 2)
}
print(arrowplot_residency_ggobj)
dev.off()

# Volcanoplot

# degs_df$degs %>%
#   data.frame() %>%
#   dplyr::mutate(label = ifelse(gene %in% genes_of_interest, gene, NA),
#                 Significance = ifelse(padj<0.05, "Significant", "NS")) %>%
#   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
#   geom_hline(yintercept = 0, col = "#BEBEBE") +
#   geom_vline(xintercept = 0, col = "#BEBEBE") +
#   geom_point_rast(aes(col = Significance)) +
#   geom_point_rast(data = . %>%
#                     dplyr::filter(!is.na(label)),
#                   col = "red") +
#   geom_label_repel(aes(label = label)) +
#   labs(title = "M1-like vs M2-like",
#        y = bquote('-'~log[10]~'(p-value)'),
#        x = bquote('-'~log[2]~'(fold change)')) +
#   scale_color_manual(values = c("NS" = "#d3d3d3", "Significant" = "#000000")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         legend.pos = "bottom")

sessionInfo()